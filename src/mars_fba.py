"""
Mars Yüzey Koşulları İçin Minimal Sentetik Hücre (JCVI-syn3A) Modellemesi
mars-minimal-cell-network — mars-minimal-gene-network'ün devamı/kardeş projesi

Bu script:
  1) JCVI-syn3A'nın yayınlanmış genom-ölçekli metabolik modelini (iMMSYN,
     155 gen, 338 reaksiyon, 304 metabolit) yükler
  2) Referans (yayınlanmış, zengin/tanımsız ortam) büyüme oranını hesaplar
  3) Mars yüzey koşullarını (O2, CO2, organik karbon, su, radyasyon->bakım
     enerjisi) sayısal kısıtlara çevirip modele uygular
  4) Mars koşulunda büyüme oranını yeniden hesaplar ve iki sonucu karşılaştırır

Model kaynağı: Breuer ve ark. 2019, eLife, "Essential metabolism for a
minimal cell" (DOI: 10.7554/eLife.36842, PMC6609329). Model dosyası
makalenin Supplementary file 9'u (SBML/FBC formatı):
https://cdn.elifesciences.org/articles/36842/elife-36842-supp9-v3.xml.zip
Model id="MMSYN", 155 gen / 338 reaksiyon / 304 metabolit — makalenin
bildirdiği sayılarla birebir doğrulandı.

ÖNEMLİ FARK (mars-minimal-gene-network'teki iYO844/iMB631'e göre):
  - Referans ortam: iYO844 için BiGG'nin tanımlı minimal ortamı kullanılmıştı.
    JCVI-syn3A için literatürde "normal büyümeyi destekleyen tanımlı bir
    ortam henüz elde edilmedi" (Breuer ve ark. 2019) -- model, yayınlanmış
    haliyle zengin/tanımsız bir ortam varsayıyor (78 exchange bileşeninin
    çoğu sınırsız alım). Bu yüzden "Dünya benzeri" referans burada modelin
    kendi varsayılan (yayınlanmış) ortamı -- bu, makalenin kendi metodolojisi
    ile tutarlı, ama iYO844/iMB631'deki tanımlı-ortam referanslarından farklı
    bir varsayım niteliğinde olduğu README'de açıkça belirtilmeli.
  - Bakım enerjisi (NGAM): iYO844'teki gibi ayrı, sabit (lb=ub) bir "ATPM"
    reaksiyonu YOK. Bunun yerine makale (Bölüm 'GAM/NGAM') NGAM'ı üç
    bileşene dağıtıyor: (1) "ATPase" reaksiyonu üzerinde tek yönlü bir alt
    sınır (0.57 mmol/gDW/h -- modelde tam olarak 0.575), (2) Protein_degrad
    alt sınırı (3.5e-4), (3) RNA_degrad alt sınırı (7.7e-3). Bu üç değer
    modelde makaledeki sayılarla birebir doğrulandı. Bu script, önceki
    projelerdeki "bakım çarpanı" yöntemini ATPase reaksiyonunun sabit
    akışına uygulayarak genişletiyor (bkz. mars_kisitlarini_uygula).
    NOT: ATPase reaksiyonu modelde TERS yönde yazılmış
    (adp+4h_e+pi <-- atp+h2o+3h_c), yani fiziksel akış negatif yönde --
    sabitleyen kısıt lower_bound değil UPPER_bound'dur (-0.575). Çarpanı
    upper_bound'a uyguluyoruz.

Windows Unicode kullanıcı adı sorunu: mars-minimal-gene-network'te
keşfedilen aynı çözüm burada da geçerli -- .xml.gz dosyasını Python'ın
gzip modülüyle açıp SBML içeriğini libSBML'e dosya yolu değil STRING
olarak veriyoruz.
"""

import gzip
import io
import os

import cobra

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_ONBELLEK = os.path.join(PROJE_KOKU, "data", "models", "iMMSYN.xml.gz")

# FBA'da çok düşük büyüme oranlarında solver'ın varsayılan feasibility
# tolerance'ı (1e-7) sahte-feasible sonuç verebilir (bkz.
# mars-minimal-gene-network/src/mars_gen_silme.py'deki kritik ders).
# Bu projede de her analizde tolerance baştan 1e-9'a çekilecek.
SOLVER_TOLERANCE = 1e-9


def modeli_yukle():
    """iMMSYN (JCVI-syn3A) modelini yerel önbellekten yükler.

    Not: dosya yolunu doğrudan cobra.io.read_sbml_model'e vermiyoruz --
    libSBML, Windows'ta yolda ASCII-dışı karakter (ör. "Ergün") olduğunda
    dosyayı açamıyor. Bunun yerine gzip'i Python'da açıp ham SBML metnini
    string olarak veriyoruz.
    """
    if not os.path.exists(MODEL_ONBELLEK):
        raise FileNotFoundError(
            f"Model önbellekte bulunamadı: {MODEL_ONBELLEK}\n"
            "data/models/iMMSYN.xml.gz eksik. Kaynak: "
            "https://cdn.elifesciences.org/articles/36842/elife-36842-supp9-v3.xml.zip "
            "(Breuer ve ark. 2019, eLife 36842, Supplementary file 9)"
        )
    with gzip.open(MODEL_ONBELLEK, "rt", encoding="utf-8") as f:
        sbml_metni = f.read()
    model = cobra.io.read_sbml_model(io.StringIO(sbml_metni))
    model.solver.configuration.tolerances.feasibility = SOLVER_TOLERANCE
    print(f"Model yüklendi: {len(model.reactions)} reaksiyon, {len(model.genes)} gen, "
          f"{len(model.metabolites)} metabolit")
    return model


def referans_buyume(model):
    """Modelin yayınlanmış (zengin/tanımsız) varsayılan ortamındaki büyüme oranı."""
    baseline = model.optimize()
    print(f"Referans (yayınlanmış ortam) büyüme oranı (1/saat): {baseline.objective_value}")
    return baseline


def bakim_reaksiyonunu_bul(model):
    """NGAM'ın ATPase bileşenini bulur ve standart yöne çevirir.

    DÜZELTME (bu projede canlı yakalanan bir gen-silme hatası): modelin
    yayınlanmış hali bu reaksiyonu TERS yönde yazıyor (adp+4h_e+pi <--
    atp+h2o+3h_c), sabitleyen kısıt upper_bound=-0.575 (zorunlu bir NEGATİF
    akış). cobra'nın standart gen-silme mekanizması bir geni sildiğinde
    ilgili reaksiyonun sınırlarını (0,0)'a çeker -- bu reaksiyon için bu,
    "kapasiteyi kaldırmak" değil "zorunlu ATP-tüketim yükünü kaldırmak"
    anlamına geliyor. Sonuç: ATPase'in tek genlerinden biri "silindiğinde"
    büyüme DÜŞMÜYOR, tam tersi ARTIYOR (0.342 -> 0.346) -- gen-silme
    analizini anlamsızlaştıran bir bulgu.
    Çözüm: reaksiyonun stokiyometrisini -1 ile çarpıp standart yöne
    çeviriyoruz (pozitif akış = fizyolojik ATP hidrolizi/bakım), sınırları
    da buna göre çeviriyoruz: eski (-1000, -0.575) -> yeni (0.575, 1000).
    Bu, artık iYO844'teki ATPM ile birebir aynı sözleşmeyi (lower_bound =
    zorunlu minimum bakım akışı) taşıyor -- gen-silme (0,0) sınırı artık
    doğru şekilde "kapasiteyi kaldırmak" anlamına geliyor.
    """
    atpase = model.reactions.get_by_id("ATPase")
    eski_sinirlar = atpase.bounds
    eski_rxn = atpase.reaction
    # Stokiyometriyi tersine çevir: her katsayıya -2x ekleyerek net -1x'e getir.
    atpase.add_metabolites({met: -2 * katsayi for met, katsayi in atpase.metabolites.items()})
    yeni_alt = -eski_sinirlar[1]
    yeni_ust = -eski_sinirlar[0]
    atpase.bounds = (yeni_alt, yeni_ust)
    print(f"Bakım (NGAM/ATPase) reaksiyonu: {atpase.id} | standart yöne çevrildi "
          f"({eski_rxn} [{eski_sinirlar}] -> {atpase.reaction} [{atpase.bounds}])")
    return atpase


def mars_kisitlarini_uygula(
    model, atpase, o2_lb=-0.5, glc_lb=-0.05, h2o_cap=1.0, bakim_carpani=3, bakim_taban=None, sessiz=False
):
    # O2: Mars atmosferinin sadece ~%0.13'ü O2 + toplam basınç Dünya'nın ~%0.6'sı -> ciddi kısıtla
    model.reactions.EX_o2_e.lower_bound = o2_lb

    # CO2: Mars atmosferinin ~%95.54'ü CO2 -> alıma açıkça izin ver, kısıtlayıcı olmasın
    model.reactions.EX_co2_e.bounds = (-1000, 1000)

    # Organik karbon: Mars yüzeyinde serbest glikoz/organik karbon yok denecek kadar az.
    # Not: modelin yayınlanmış hali EX_glc__D_e için zaten bir üst sınır (-7.4) içeriyor
    # (muhtemelen proteomik-tabanlı taşıyıcı kapasitesi); Mars kısıtı bunu daha da daraltır.
    model.reactions.EX_glc__D_e.lower_bound = glc_lb

    # Su: düşük su aktivitesi (aw ~0.4) -> su akışını daralt (hem alım hem metabolik atım).
    model.reactions.EX_h2o_e.bounds = (-h2o_cap, h2o_cap)

    # Bakım enerjisi: radyasyon hasarını onarmak ek ATP gerektirir -> NGAM'ın ATPase
    # bileşenini "bakim_carpani" kat artırıyoruz. bakim_reaksiyonunu_bul() artık
    # ATPase'i standart yöne çevirdiği için (bkz. o fonksiyonun docstring'i),
    # sabitleyen kısıt artık iYO844'teki ATPM ile birebir aynı: LOWER_bound
    # (pozitif, zorunlu minimum akış).
    #
    # DİKKAT (mars-minimal-gene-network'te öğrenilen bir üretim hatasından ders):
    # çarpanı reaksiyonun O ANKİ sınırına uygularsak, aynı model/atpase nesnesi
    # birden çok kez (ör. bir tarama döngüsünde) bu fonksiyondan geçirildiğinde
    # değer her seferinde katlanarak büyür. Bunu önlemek için taban değeri HER
    # ZAMAN açıkça bilinen bir referanstan alınır: çağıran "bakim_taban" vermezse,
    # atpase.lower_bound SADECE bu fonksiyon hiç çağrılmamışsa (yani hâlâ modelin
    # orijinal değeriyse) taban olarak kabul edilir -- modeli tekrar kullanan
    # çağıranlar (ör. duyarlılık/gen silme analizi) bakim_taban'ı MUTLAKA açıkça
    # vermeli.
    taban = bakim_taban if bakim_taban is not None else atpase.lower_bound
    yeni_alt_sinir = taban * bakim_carpani
    atpase.bounds = (yeni_alt_sinir, 1000)
    if not sessiz:
        print(f"Yeni bakım (NGAM/ATPase) sınırı: {atpase.bounds}")

    return model


def mars_buyume(model):
    """Mars koşulundaki FBA çözümünü döndürür.

    DİKKAT (bu projede canlı yakalanmış bir solver-artefaktı): model
    infeasible olduğunda glpk/optlang bazen sıfır olmayan, anlamsız bir
    objective_value döndürebiliyor -- özellikle aynı model nesnesi üzerinde
    daha önce (ör. referans_buyume ile) bir optimize() çağrısı yapılmışsa,
    warm-start solver durumu bir sonraki (infeasible) çağrıya taşınabiliyor.
    Bu yüzden status EXPLICIT olarak kontrol ediliyor; infeasible ise
    objective_value hiçbir zaman gerçek bir büyüme oranı olarak
    raporlanmıyor/kullanılmıyor.
    """
    mars_solution = model.optimize(raise_error=False)
    durum = model.solver.status
    if durum != "optimal":
        print(f"Mars koşulunda büyüme oranı: TANIMSIZ (durum: {durum} -- model infeasible, "
              "objective_value güvenilir değil, raporlanmıyor)")
    else:
        print(f"Mars koşulunda büyüme oranı (1/saat): {mars_solution.objective_value} | durum: {durum}")
    return mars_solution


def main():
    model = modeli_yukle()
    baseline = referans_buyume(model)
    atpase = bakim_reaksiyonunu_bul(model)
    model = mars_kisitlarini_uygula(model, atpase)
    mars_solution = mars_buyume(model)

    print()
    print("--- Özet ---")
    print(f"Referans (yayınlanmış ortam): {baseline.objective_value}")
    if model.solver.status == "optimal":
        print(f"Mars koşulu:                  {mars_solution.objective_value}")
    else:
        print(f"Mars koşulu:                  TANIMSIZ (infeasible)")


if __name__ == "__main__":
    main()
