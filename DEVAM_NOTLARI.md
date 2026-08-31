# Proje Devam Notları — mars-minimal-cell-network

Son güncelleme: 2026-08-31

Bu dosya `mars-minimal-gene-network/DEVAM_NOTLARI.md`'nin devamıdır — o
dosyanın madde 5-6'sında planlanan "minimal sentetik hücre" yönü burada
ayrı bir proje olarak yürütülüyor. Ana repo ve tüm geçmiş bulgular için
önce oraya bak (özellikle Windows Unicode-path sorunu, solver tolerance
dersi, git commit/push blanket onayı gibi genel dersler hâlâ geçerli).

## Şu ana kadar yapılanlar

1. **GEM taraması bitirildi**: JCVI-syn3A için gerçek, yayınlanmış bir
   metabolik model bulundu — Breuer ve ark. 2019, eLife 36842,
   Supplementary file 9 (SBML). `Luthey-Schulten-Lab/Minimal_Cell`
   reposu incelendi ama o whole-cell kinetik simülasyon için (FBA modeli
   orada değiştirilmiş/gömülü halde) — ham SBML kaynağı olarak eLife
   supplementary'si tercih edildi. iJL208 (Mesoplasma florum, 208 gen)
   farklı bir organizma olduğu için elenmedi ama kullanılmadı.
2. **Model indirildi, doğrulandı**: `data/models/iMMSYN.xml.gz`.
   155 gen / 338 reaksiyon / 304 metabolit — makaleyle birebir. Varsayılan
   ortamda büyüme 0.342/saat — makalenin ~2 saat çiftlenme süresiyle
   örtüşüyor.
3. **NGAM yapısı literatürle doğrulandı**: makalenin tam metni (JATS XML,
   `cdn.elifesciences.org/articles/36842/elife-36842-v3.xml`) çekilip
   'GAM/NGAM' bölümü okundu. Modeldeki `ATPase` (0.575), `Protein_degrad`
   (0.00035), `RNA_degrad` (0.0077) sabit alt sınırları makaledeki
   sayılarla birebir eşleşti. Ad hoc bir bakım-enerjisi varsayımı
   uydurmaya gerek kalmadı.
4. **`src/mars_fba.py` yazıldı** — B. subtilis projesindeki yapıyı
   (modeli_yukle/referans_buyume/bakim_reaksiyonunu_bul/
   mars_kisitlarini_uygula/mars_buyume) aynen taşıyor, NGAM/ATPase'in
   ters yön sözleşmesine göre uyarlanmış.
5. **Solver artefaktı canlı yakalandı ve düzeltildi**: `main()` içinde
   art arda `optimize()` çağrıları, infeasible bir çözüm için glpk'nin
   sıfır olmayan anlamsız bir `objective_value` döndürmesine yol açtı
   (warm-start durumu taşınması). `mars_buyume()` artık
   `model.solver.status`'u EXPLICIT kontrol ediyor, infeasible ise
   objective_value hiç raporlanmıyor. **Bu türden bir kontrolü her yeni
   FBA script'inde tekrarla.**
6. **İlk (duyarlılık analizinden önce) bulgu**: B. subtilis parametreleri
   birebir taşındığında (o2=-0.5, glc=-0.05, h2o=1.0, çarpan=3) model
   infeasible. Kısıtlar izole edilerek tarandığında: **su B.
   subtilis'teki gibi tek başına kısıtlayıcı DEĞİL** burada; asıl
   kısıtlayıcı **glikoz** — `EX_glc__D_e` alt sınırında -0.8 (feasible,
   büyüme≈0.0016/saat) ile -0.75 (infeasible) arasında keskin bir uçurum
   var. Detay ve yorum için README > "Şu ana kadarki bulgu".

## Henüz yapılmadı / sıradaki somut adımlar

1. **Tam duyarlılık analizi** — B. subtilis/Salinibacter projelerindeki
   `mars_duyarlilik.py` şablonunu al, dört parametreyi (o2_lb, glc_lb,
   h2o_cap, bakim_carpani) tarayan bir `src/mars_duyarlilik.py` yaz.
   Özellikle glc_lb ekseni etrafında ince taramaya odaklan (uçurum -0.8/
   -0.75 civarında).
2. **Gen esansiyellik/silme analizi** — `src/mars_gen_silme.py` şablonunu
   al, `SOLVER_TOLERANCE = 1e-9`'u BAŞTAN uygula (bkz. ana projedeki
   kritik ders). Feasible bir Mars parametre noktası seçmek gerekecek
   (ör. glc_lb=-0.8 civarı, infeasible noktada essentiality analizi
   anlamsız).
3. **Üç modelin karşılaştırılması**: B. subtilis (su-kısıtlı uçurum),
   Salinibacter (uçurum yok, doğrusal), JCVI-syn3A (glikoz-kısıtlı
   uçurum) — üç farklı "kısıtlayıcı darboğaz" profili ortaya çıktı. Bu,
   makalenin ana karşılaştırma bulgusu olabilir; ayrı bir
   `karsilastirma.py`/tablo ile bir araya getirilmeli.
4. **GitHub'a push** — yerel git init/commit yapıldıktan sonra kullanıcıya
   yeni bir GitHub reposu (`mars-minimal-cell-network`) oluşturup push
   etme izni ayrıca soruldu/soruluyor (önceki blanket onay mevcut bir
   repoya commit/push için verilmişti, yeni bir public repo açmak farklı
   bir eylem olarak ayrıca teyit edildi).
5. **README'de dürüstlük notu eklenmeli**: referans ortamın (zengin/
   tanımsız) B. subtilis/Salinibacter'deki tanımlı-ortam referanslarından
   niteliksel farkı, makalede/README'de MUTLAKA vurgulanmalı — aksi halde
   üç model arası büyüme oranı karşılaştırması yanıltıcı olur.

## Genel hatırlatmalar (ana projeden taşınan, hâlâ geçerli)

- Windows Unicode kullanıcı adı → gzip+string SBML yükleme yöntemi.
- Solver tolerance: essentiality/gen silme analizlerinde MUTLAKA
  `1e-9`'a çek.
- FBA sonucunu asla `objective_value`'ya körü körüne güvenerek raporlama
  — `model.solver.status`'u her zaman kontrol et (bu projede madde 5'te
  tekrar doğrulandı).
- Kullanıcı Türkçe konuşuyor, dürüst/kaynaklı/"bulunamadı"yı da rapor
  eden üslup — literatür atıflarını asla hafızadan yazma, her zaman
  WebSearch/WebFetch ile çapraz doğrula.
