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

from mars_fba import modeli_yukle
from mars_gen_silme import NGAM_ISTISNA_GENLERI

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


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    print("=== JCVI-syn3A (iMMSYN) -- kısıtsız/yayınlanmış referans için ARDIŞIK indirgeme ===")
    model = modeli_yukle()
    tutulan, cikarilan, son_sol = ardisik_indirgeme(model, korunacak_genler=NGAM_ISTISNA_GENLERI)

    print(f"\nToplam gen: {len(tutulan) + len(cikarilan)}")
    print(f"Tutulan (gerçekten gerekli) gen: {len(tutulan)}")
    print(f"Çıkarılan (gerçekten gereksiz) gen: {len(cikarilan)}")
    print(f"İndirgenmiş ağ son durum: {model.solver.status}, büyüme: {son_sol.objective_value}")

    eski_yol = os.path.join(SONUC_KLASORU, "gen_silme_sonuclari.csv")
    if os.path.exists(eski_yol):
        eski = pd.read_csv(eski_yol)
        eski_esansiyel = set(eski[(eski.senaryo == "Referans_kisitsiz") & (eski.esansiyel == True)].gen_id)
        print(f"\nKarşılaştırma: tekli-silme (düzeltilmiş) esansiyel sayısı: {len(eski_esansiyel)}")
        print(f"Bu script'in tuttuğu gen sayısı: {len(tutulan)}")

    pd.DataFrame({"gen_id": tutulan, "durum": "tutuldu_gerekli"}).to_csv(
        os.path.join(SONUC_KLASORU, "minimal_ag_tutulan_genler.csv"), index=False)
    pd.DataFrame({"gen_id": cikarilan, "durum": "cikarildi_gereksiz"}).to_csv(
        os.path.join(SONUC_KLASORU, "minimal_ag_cikarilan_genler.csv"), index=False)
    print("\nKaydedildi: results/minimal_ag_*.csv")


if __name__ == "__main__":
    main()
