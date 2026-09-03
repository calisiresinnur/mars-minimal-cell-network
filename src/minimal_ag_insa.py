"""
GERÇEK minimal gen ağı inşası — JCVI-syn3A (iMMSYN), ardışık indirgeme.

Bkz. mars-minimal-gene-network/src/minimal_ag_insa.py (B. subtilis
versiyonu, aynı yöntem/gerekçe) — kullanıcının kritik sorusu üzerine:
tekli gen silme "esansiyel" listesi izoenzim/yedek-yol gruplarını
yakalayamıyor. Bu script genleri TEK TEK, o ana kadar indirgenmiş ağa
göre test ederek gerçekten işlevsel bir minimal ağ inşa ediyor.

ÖZEL DURUM (bu proje için) -- NGAM istisna genleri: `mars_gen_silme.py`
9 geni (ATPase×8 + Protein_degrad×1) "silinince büyüme ARTIYOR" paradoksu
nedeniyle ELLE esansiyel olarak düzeltmişti (makale Table 4 ile çapraz
doğrulanarak). Bu greedy indirgeme algoritması, düzeltme uygulanmazsa
AYNI paradoksa düşüp bu 9 geni de "gereksiz" diye çıkarır -- bu yüzden
bu 9 gen BAŞTAN "tutulan" listesine ekleniyor, hiç test edilmiyor.
"""

import os

import pandas as pd

from mars_fba import bakim_reaksiyonunu_bul, mars_kisitlarini_uygula, modeli_yukle
from mars_gen_silme import MARS_SENARYOLARI, NGAM_ISTISNA_GENLERI

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SONUC_KLASORU = os.path.join(PROJE_KOKU, "results")

ESANSIYELLIK_ESIGI = 0.01
SOLVER_TOLERANCE = 1e-9


def ardisik_indirgeme(model, esik=ESANSIYELLIK_ESIGI, korunacak_genler=None, sessiz=False):
    model.solver.configuration.tolerances.feasibility = SOLVER_TOLERANCE
    wt = model.optimize(raise_error=False)
    if model.solver.status != "optimal":
        raise RuntimeError(f"WT durumu optimal değil ({model.solver.status})")
    wt_buyume = wt.objective_value
    if not sessiz:
        print(f"WT büyüme: {wt_buyume:.6f}")

    korunacak_genler = set(korunacak_genler or [])
    tum_gen_id = [g.id for g in model.genes]
    cikarilan, tutulan = [], []

    for i, gid in enumerate(tum_gen_id):
        gene = model.genes.get_by_id(gid)
        if not gene.functional:
            continue
        if gid in korunacak_genler:
            # bkz. modül docstring'i -- NGAM istisna genleri hiç test edilmiyor,
            # doğrudan tutuluyor (test edilirse paradoks nedeniyle yanlışlıkla
            # "gereksiz" çıkarlardı).
            tutulan.append(gid)
            continue
        with model:
            gene.knock_out()
            sol = model.optimize(raise_error=False)
            buyume = sol.objective_value if model.solver.status == "optimal" else 0.0
        oran = (buyume / wt_buyume) if wt_buyume else 0.0
        if oran >= esik:
            gene.knock_out()
            cikarilan.append(gid)
        else:
            tutulan.append(gid)
        if not sessiz and (i + 1) % 50 == 0:
            print(f"  {i + 1}/{len(tum_gen_id)} -- tutulan: {len(tutulan)}, çıkarılan: {len(cikarilan)}")

    son_sol = model.optimize(raise_error=False)
    return tutulan, cikarilan, son_sol


def senaryo_calistir(etiket, kisit_uygula):
    model = modeli_yukle()
    if kisit_uygula is not None:
        kisit_uygula(model)
    print(f"\n=== {etiket} ===")
    tutulan, cikarilan, son_sol = ardisik_indirgeme(model, korunacak_genler=NGAM_ISTISNA_GENLERI)
    print(f"Toplam gen: {len(tutulan) + len(cikarilan)} | Tutulan: {len(tutulan)} | "
          f"Çıkarılan: {len(cikarilan)} | Son durum: {model.solver.status} | "
          f"Büyüme: {son_sol.objective_value}")
    return etiket, set(tutulan), set(cikarilan), son_sol.objective_value if model.solver.status == "optimal" else None


def mars_kisiti(senaryo):
    def uygula(model):
        atpm = bakim_reaksiyonunu_bul(model)
        mars_kisitlarini_uygula(
            model, atpm, o2_lb=senaryo["o2"], glc_lb=senaryo["glc"], h2o_cap=senaryo["h2o"],
            bakim_carpani=senaryo["bakim_carpani"], sessiz=True,
        )
    return uygula


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    sonuclar = [senaryo_calistir("Referans_kisitsiz", None)]
    for s in MARS_SENARYOLARI:
        sonuclar.append(senaryo_calistir(s["etiket"], mars_kisiti(s)))

    eski = pd.read_csv(os.path.join(SONUC_KLASORU, "gen_silme_sonuclari.csv"))

    print("\n\n=== ÖZET: her senaryo için tekli-silme vs gerçek minimal ağ ===")
    satirlar = []
    for etiket, tutulan, cikarilan, buyume in sonuclar:
        eski_esansiyel = set(eski[(eski.senaryo == etiket) & (eski.esansiyel == True)].gen_id)
        print(f"{etiket:24s}: tekli-silme={len(eski_esansiyel):4d}  gerçek_minimal={len(tutulan):4d}  "
              f"fark=+{len(tutulan - eski_esansiyel):3d}  büyüme={buyume}")
        satirlar.append(dict(senaryo=etiket, tekli_silme_sayisi=len(eski_esansiyel),
                              gercek_minimal_sayisi=len(tutulan), fark=len(tutulan - eski_esansiyel),
                              indirgenmis_ag_buyume=buyume))
        pd.DataFrame({"gen_id": sorted(tutulan), "durum": "tutuldu_gerekli"}).to_csv(
            os.path.join(SONUC_KLASORU, f"minimal_ag_tutulan_genler_{etiket}.csv"), index=False)

    pd.DataFrame(satirlar).to_csv(os.path.join(SONUC_KLASORU, "minimal_ag_ozet_tum_senaryolar.csv"), index=False)

    print("\n=== Referans minimal ağı ile Mars minimal ağları arasındaki fark ===")
    ref_tutulan = sonuclar[0][1]
    for etiket, tutulan, _, _ in sonuclar[1:]:
        sadece_mars = tutulan - ref_tutulan
        sadece_ref = ref_tutulan - tutulan
        print(f"{etiket}: Sadece Mars'ta gerekli: {len(sadece_mars)} gen "
              f"({', '.join(sorted(sadece_mars)) if sadece_mars else '-'}) | "
              f"Sadece referansta gerekli: {len(sadece_ref)} gen")

    print("\nKaydedildi: results/minimal_ag_*.csv")


if __name__ == "__main__":
    main()
