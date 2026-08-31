"""
Tekli gen silme (single gene deletion) analizi — JCVI-syn3A (iMMSYN).

mars-minimal-gene-network/src/mars_gen_silme.py'nin yöntemini izliyor:
referans (kısıtsız/yayınlanmış ortam) + üç Mars senaryosu için modelin
155 geninin her biri tek tek "silinip" büyüme oranı yeniden hesaplanıyor.

NEDEN t*'A ÇOK YAKIN BİR MARJ DEĞİL, %5+ BÜYÜME NOKTASI (mars-minimal-
gene-network'teki t*+0.01 stratejisinden BİLEREK farklı bir karar):
mars_duyarlilik.py ile tam feasibility sınırını (t*) bisection ile
ararken CANLI bir numerik kırılganlık yakalandı -- t* değerini 6 ondalık
basamağa yuvarlayıp yeni bir script'e elle kopyaladığımda (~1e-6
hassasiyet kaybı) AYNI NOKTA "optimal" yerine "infeasible" çıktı. Daha da
çarpıcısı: sınıra çok yakın bir ilk denemede (t*+0.04, WT büyüme ~%2
referans), `single_gene_deletion` çıktısında bazı KO'lar için (gerçekte
infeasible olmaları gerekirken) daha önce mars_fba.py'de canlı yakalanan
AYNI sahte-değer (0.0013440550524949...) tekrar ortaya çıktı -- yani bu
modelde düşük-büyüme Mars senaryolarında solver artefaktları B.
subtilis'tekinden DAHA agresif. Bu yüzden burada çok daha rahat bir marj
seçildi: WT büyüme referansın en az %5'i olacak şekilde t taranmış, HER
seçilen nokta kullanılmadan önce 5 kez BAĞIMSIZ (taze model, sıfırdan
optimize) tekrarla "optimal" durumun bit-bit aynı çıktığı doğrulanmıştır.

DÜZELTME 2 -- NGAM-ilişkili genlerin essentiality'si (kullanıcının
"mantık olarak da denetle" talebiyle bulundu, makale Table 4 ile
çapraz kontrol edilerek): `ATPase`'e bağlı 8 gen (MMSYN1_0789-0796,
F1F0-ATP sentaz alt birimleri) VE `Protein_degrad`'a bağlı 1 gen
(MMSYN1_0394) -- toplam 9 gen -- standart cobra gen-silmesinde
"silindiğinde" büyüme ARTIYOR (bounds (0,0)'a çekiliyor, bu da
kapasiteyi değil zorunlu NGAM YÜKÜNÜ kaldırıyor). İlk yorumum "bu bir
model sınırlaması, gerçek biyolojiyle çelişebilir ama hesaplama hatası
değil" şeklindeydi. AMA makalenin kendi Table 4'ü (locus-düzeyinde
Ess_FBA sütunu) bu 9 genin TAMAMINI ■ (esansiyel) olarak işaretliyor --
yani makalenin KENDİ in silico essentiality metodolojisi bu genleri
esansiyel buluyor, benim naif cobra knockout'um değil.

Doğru yorum: NGAM'ın dayattığı minimum akış (0.575 ATPase, 0.00035
Protein_degrad) hücrenin GERÇEK, sürekli bir fizyolojik ihtiyacını
(pH homeostazı, protein döngüsü) temsil ediyor -- bu ihtiyaç, ilgili
enzim geni silindiğinde ORTADAN KALKMAZ, sadece KARŞILANAMAZ hale gelir.
cobra'nın standart knockout'u (bounds->(0,0)) bunu "artık bu kadar ATP
harcamana gerek yok" gibi yanlış yorumluyor; doğrusu "artık bu ihtiyacı
karşılayamıyorsun -> infeasible" olmalı. Bu yüzden bu 9 gen, ham
simülasyon sonucundan BAĞIMSIZ olarak esansiyel kabul ediliyor (bkz.
NGAM_ISTISNA_GENLERI ve main()'deki düzeltme).

DOĞRULAMA: bu düzeltme sonrası referans esansiyel gen sayısı
114 (ham) + 9 (düzeltme) = 123/155 (%79.4) -- makalenin bildirdiği
123/155 (%79) ile TAM ÖRTÜŞÜYOR. Ayrıca birkaç bağımsız gen (ALATRS/
0163, ARGTRS/0535, ASNTRS/0076, ASPTRS/0287, CYSTRS/0837 -- hepsi
tRNA sentetaz, makalede Ess_FBA=■) benim ham sonuçlarımla da zaten
örtüşüyordu -- yani genel pipeline baştan beri doğruydu, sadece bu 9
NGAM-geni özel durumu atlanmıştı.

Çıktı:
  - results/gen_silme_sonuclari.csv          (senaryo x gen x büyüme, ham veri)
  - results/mars_yeni_esansiyel_genler.csv   (Mars'ta esansiyel, referansta değil)
  - results/mars_dispanse_olan_genler.csv    (referansta esansiyel, Mars'ta değil)
"""

import os

import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from mars_fba import bakim_reaksiyonunu_bul, mars_kisitlarini_uygula, modeli_yukle

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SONUC_KLASORU = os.path.join(PROJE_KOKU, "results")

ESANSIYELLIK_ESIGI = 0.01  # KO büyüme / WT büyüme bu değerin altındaysa "esansiyel"
ISLEMCI_SAYISI = 1  # bkz. modül docstring'i -- bu modelde düşük-büyüme noktalarında
# solver artefaktları çok agresif; tek işlemci tercih edildi (155 gen için hızlı zaten,
# multiprocessing'in ek bir pickling/state riski katmasına gerek yok)
SOLVER_TOLERANCE = 1e-9  # bkz. mars_fba.py / mars-minimal-gene-network'teki kritik ders

# NGAM'a bağlı, standart knockout'ta "büyüme artıyor" paradoksu gösteren genler --
# makalenin Table 4'üyle çapraz kontrol edilerek esansiyel kabul ediliyor
# (bkz. modül docstring'indeki "DÜZELTME 2").
ATPASE_GENLERI = [
    "MMSYN1_0789", "MMSYN1_0790", "MMSYN1_0791", "MMSYN1_0792",
    "MMSYN1_0793", "MMSYN1_0794", "MMSYN1_0795", "MMSYN1_0796",
]
NGAM_ISTISNA_GENLERI = ATPASE_GENLERI + ["MMSYN1_0394"]  # + Protein_degrad geni

# mars_duyarlilik.py sonuçlarından: WT büyüme referansın en az %5'i olacak şekilde
# seçilmiş t değerleri -- her biri 5x bağımsız tekrarla "optimal" tutarlılığı
# doğrulandı (bkz. modül docstring'i). SERT/ILIMLI ankorları mars_duyarlilik.py ile aynı.
MARS_SENARYOLARI = [
    dict(etiket="Mars_bakim_x1.5_marj", t=0.54, o2=-5.6300, glc=-1.1030, h2o=16.6600, bakim_carpani=1.5),
    dict(etiket="Mars_bakim_x2.0_marj", t=0.58, o2=-6.0100, glc=-1.1810, h2o=17.8200, bakim_carpani=2.0),
    dict(etiket="Mars_bakim_x3.0_marj", t=0.66, o2=-6.7700, glc=-1.3370, h2o=20.1400, bakim_carpani=3.0),
]


def senaryo_calistir(etiket, kisit_uygula):
    """kisit_uygula(model) -> None; None geçilirse referans (kısıtsız/yayınlanmış ortam)."""
    model = modeli_yukle()  # zaten SOLVER_TOLERANCE=1e-9 ile geliyor (bkz. mars_fba.py)
    if kisit_uygula is not None:
        kisit_uygula(model)

    wt = model.optimize(raise_error=False)
    print(f"{etiket}: WT büyüme = {wt.objective_value} (durum: {model.solver.status})")
    if model.solver.status != "optimal":
        raise RuntimeError(f"{etiket}: WT durumu optimal değil ({model.solver.status}) -- senaryo geçersiz")

    sonuc = single_gene_deletion(model, processes=ISLEMCI_SAYISI)
    sonuc = sonuc.reset_index(drop=True)
    sonuc["gen_id"] = sonuc["ids"].apply(lambda s: next(iter(s)) if s else None)
    sonuc["senaryo"] = etiket
    sonuc["wt_buyume"] = wt.objective_value
    # DÜZELTME (bu projede canlı yakalanan ikinci bir hata): status='infeasible'
    # olan KO'larda cobra "growth" alanını NaN bırakıyor. "oran = growth/wt_buyume"
    # de NaN olur, ve "oran < ESANSIYELLIK_ESIGI" pandas'ta NaN için HER ZAMAN
    # False döner -- yani infeasible (en kesin esansiyel durum!) YANLIŞLIKLA
    # "esansiyel DEĞİL" sayılıyordu. İlk çalıştırmada bu, "0 esansiyel gen"
    # gibi biyolojik olarak imkânsız bir sonuca yol açtı. DÜZELTİLMİŞ mantık:
    # status != 'optimal' ise büyüme=0 kabul edilir (KO modeli tamamen
    # infeasible yapıyor -> kesin esansiyel).
    sonuc["growth_efektif"] = sonuc["growth"].where(sonuc["status"] == "optimal", 0.0)
    sonuc["oran"] = sonuc["growth_efektif"] / wt.objective_value
    sonuc["esansiyel"] = sonuc["oran"] < ESANSIYELLIK_ESIGI
    return sonuc[["gen_id", "senaryo", "growth", "wt_buyume", "oran", "status", "esansiyel"]]


def mars_kisiti(senaryo):
    def uygula(model):
        atpase = bakim_reaksiyonunu_bul(model)
        mars_kisitlarini_uygula(
            model, atpase,
            o2_lb=senaryo["o2"], glc_lb=senaryo["glc"], h2o_cap=senaryo["h2o"],
            bakim_carpani=senaryo["bakim_carpani"], sessiz=True,
        )
    return uygula


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    tum_sonuclar = [senaryo_calistir("Referans_kisitsiz", None)]
    for s in MARS_SENARYOLARI:
        tum_sonuclar.append(senaryo_calistir(s["etiket"], mars_kisiti(s)))

    df = pd.concat(tum_sonuclar, ignore_index=True)
    df["atpase_geni_mi"] = df["gen_id"].isin(ATPASE_GENLERI)
    df["ngam_istisnasi_mi"] = df["gen_id"].isin(NGAM_ISTISNA_GENLERI)
    # DÜZELTME 2 (bkz. modül docstring'i): NGAM-ilişkili genler için ham simülasyon
    # sonucu (esansiyel_ham) korunuyor ama asıl kullanılan "esansiyel" sütunu, makale
    # Table 4 ile doğrulanmış düzeltmeyi yansıtıyor.
    df["esansiyel_ham"] = df["esansiyel"]
    df["esansiyel"] = df["esansiyel_ham"] | df["ngam_istisnasi_mi"]
    csv_yolu = os.path.join(SONUC_KLASORU, "gen_silme_sonuclari.csv")
    df.to_csv(csv_yolu, index=False)
    print(f"\nHam sonuçlar kaydedildi: {csv_yolu} ({len(df)} satır)")

    print("\n--- NGAM istisna genleri (ham sonuç vs düzeltilmiş, bkz. modül docstring'i) ---")
    print(df[df["ngam_istisnasi_mi"]][["gen_id", "senaryo", "growth", "wt_buyume", "oran", "esansiyel_ham", "esansiyel"]].to_string(index=False))

    pivot = df.pivot(index="gen_id", columns="senaryo", values="esansiyel")
    referans_esansiyel = pivot["Referans_kisitsiz"]
    mars_kolonlari = [c for c in pivot.columns if c != "Referans_kisitsiz"]

    mars_yeni_esansiyel = pivot[(~referans_esansiyel) & pivot[mars_kolonlari].any(axis=1)].copy()
    mars_yeni_esansiyel["kac_mars_senaryosunda"] = mars_yeni_esansiyel[mars_kolonlari].sum(axis=1)
    mars_yeni_esansiyel = mars_yeni_esansiyel.sort_values("kac_mars_senaryosunda", ascending=False)
    mars_yeni_esansiyel.to_csv(os.path.join(SONUC_KLASORU, "mars_yeni_esansiyel_genler.csv"))

    mars_dispanse = pivot[referans_esansiyel & (~pivot[mars_kolonlari]).any(axis=1)].copy()
    mars_dispanse["kac_mars_senaryosunda_dispanse"] = (~mars_dispanse[mars_kolonlari]).sum(axis=1)
    mars_dispanse = mars_dispanse.sort_values("kac_mars_senaryosunda_dispanse", ascending=False)
    mars_dispanse.to_csv(os.path.join(SONUC_KLASORU, "mars_dispanse_olan_genler.csv"))

    print("\n--- Özet ---")
    for kolon in pivot.columns:
        print(f"{kolon:32s}: {int(pivot[kolon].sum()):4d} esansiyel gen / {len(pivot)}")
    print(f"\nMars'a özgü YENİ esansiyel gen adayı: {len(mars_yeni_esansiyel)}")
    print(f"Mars'ta esansiyel OLMAKTAN ÇIKAN gen: {len(mars_dispanse)}")
    if len(mars_yeni_esansiyel) > 0:
        print("  Yeni esansiyel:", ", ".join(mars_yeni_esansiyel.index.tolist()))
    if len(mars_dispanse) > 0:
        print("  Dispanse olan:", ", ".join(mars_dispanse.index.tolist()))
    print("Kaydedildi: results/mars_yeni_esansiyel_genler.csv, results/mars_dispanse_olan_genler.csv")


if __name__ == "__main__":
    main()
