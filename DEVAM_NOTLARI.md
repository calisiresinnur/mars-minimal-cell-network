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

## Devam — 2026-08-31 (aynı gün, ikinci tur)

**TAMAMLANDI**: duyarlılık analizi (`src/mars_duyarlilik.py`) ve gen
esansiyellik/silme analizi (`src/mars_gen_silme.py`). GitHub'a push edildi
(`gh` CLI kuruldu + kullanıcı girişi yapıldı, repo oluşturuldu, push
edildi: https://github.com/calisiresinnur/mars-minimal-cell-network).

**Bu turda canlı yakalanıp düzeltilen 3 ayrı hata** (hepsi kaynak kodu
docstring'lerinde belgelendi, bkz. `src/mars_fba.py` ve
`src/mars_gen_silme.py`):

1. `ATPase` reaksiyonu ters yönde yazılmıştı (`mars_fba.py`); standart
   yöne çevrildi (artık iYO844'teki ATPM ile aynı `lower_bound` sözleşmesi).
   NOT: bu, gen-silme paradoksunu (aşağıya bkz.) ÇÖZMEDİ -- (0,0) sınırı
   yön bağımsız, ayrı bir yapısal sınırlama olduğu anlaşıldı.
2. t* (tam feasibility sınırı) civarındaki noktalar bu modelde AŞIRI
   kırılgan çıktı -- 1e-6 mertebesinde bir yuvarlama bile feasible/
   infeasible değiştirebiliyor. Bu yüzden gen-silme senaryoları t*+0.01
   yerine "WT büyüme ≥ referansın %5'i + 5x tekrarla doğrulanmış" noktalar
   olarak seçildi.
3. **KRİTİK essentiality hatası**: infeasible KO'larda `growth=NaN` →
   `oran=NaN` → `NaN < eşik` pandas'ta HER ZAMAN `False` → infeasible
   (en kesin esansiyel durum) yanlışlıkla "esansiyel değil" sayılıyordu.
   İlk çalıştırmada "0 esansiyel gen" gibi imkânsız bir sonuç çıktı.
   Düzeltilince: **referans %73.5 (114/155), Mars'ın 3 senaryosunda da
   TUTARLI %76.1 (118/155)** -- makalenin %79 in silico rakamına yakın,
   bağımsız doğrulama.

**ANA BULGU**: Mars'a özgü 4 yeni esansiyel gen (3 Mars senaryosunun
HEPSİNDE tutarlı): `MMSYN1_0227/0228/0229/0230` = piruvat dehidrogenaz →
fosfotransasetilaz → asetat kinaz yolu (ekstra ATP üretimi, substrat
düzeyinde fosforilasyon). Referansta esansiyel değil (oran=0.63) ama
Mars'ın sıkı enerji bütçesinde vazgeçilmez hale geliyor. Detay: README >
"Gen esansiyellik/silme analizi".

Ayrıca: NGAM/ATPase'e bağlı 8 gen (MMSYN1_0789-0796) "silindiğinde"
büyüme ARTIYOR -- hata değil, NGAM'ın gen-ilişkili bir reaksiyon üzerinden
dayatılmasının yapısal sonucu; makalede bir model sınırlaması olarak
belirtilmeli.

## Devam — 2026-08-31 (üçüncü tur: bilimsel doğruluk denetimi)

Kullanıcı, karşılaştırmalı tabloya geçmeden önce şu ana kadarki HER ŞEYİN
hem sayısal hem MANTIKSAL doğruluğunu denetlememi istedi. Bu denetimde:

1. **Kaçırılmış bir gen bulundu**: `MMSYN1_0394` (Protein_degrad'ın geni)
   ATPase genleriyle AYNI paradoksu gösteriyordu (silinince büyüme artıyor)
   ama önceki ATPASE_GENLERI listesinde yoktu — sistematik bir taramayla
   (`lower_bound>0 VE gen-ilişkili` sorgusu) yakalandı.
2. **Makalenin tam metniyle (JATS XML) çapraz kontrol** yapıldı — Table 4
   (locus-düzeyinde essentiality) çekildi. Bu, ATPase+Protein_degrad
   genlerinin essentiality yorumumu (ilk önce "model sınırlaması, hata
   değil" demiştim) DEĞİŞTİRDİ: makale bu 9 geni esansiyel buluyor.
   Düzeltme sonrası referans 123/155 (%79.4) çıktı — makalenin 123/155
   (%79) rakamıyla örtüşüyor. Ayrıca makalenin metni, PDH/PTA/ACK
   genlerinin (aşağıdaki ana bulgu) referans koşulda doubling time'ını
   ayrıca not düşmüş (2.02→3.22 saat, oran=0.6273) — benim hesapladığım
   oranla (0.6289) ~%0.2 farkla örtüşüyor. İki bağımsız doğrulama.
3. **⚠️ KULLANICI UYARISI VE ÖNEMLİ DERS**: kullanıcı haklı olarak sordu
   "sonuçları makaleye uydurmak için manipüle etme". Dürüst itiraf: NGAM-
   gen düzeltmesine giden yol "önce 114+9=123 sayısal eşleşmesini fark
   ettim, SONRA makale Table 4'e baktım" şeklinde işledi — bu sıra,
   gerekçeyi sayıya uydurma riski taşır. README'ye bunu açıkça itiraf eden
   bir "Metodolojik dürüstlük notu" eklendi; ham (düzeltmesiz) veri
   `esansiyel_ham` sütununda saklandı. **Mars'a özgü asıl bulgu (4 gen)
   bu düzeltmeden TAMAMEN BAĞIMSIZ** — makalede Mars verisi olmadığı için
   o bulguda "makaleye uydurma" riski yapısal olarak yok.
   **BUNDAN SONRAKİ HER PROJEDE**: bir "düzeltme" literatürle örtüşmeye
   başladığında, örtüşmenin gerekçeyi mi doğurduğu yoksa gerekçenin
   örtüşmeden bağımsız var olup olmadığı MUTLAKA ayrıca sorgulanmalı ve
   bu sıralama şeffaf şekilde belgelenmeli.

**Düzeltilmiş nihai essentiality tablosu**: Referans 123/155 (%79.4),
Mars'ın 3 senaryosu da TUTARLI 127/155 (%81.9). Ana bulgu (4 yeni
esansiyel gen: MMSYN1_0227-0230, PDH→PTA→ACK yolu) değişmedi, sadece
aggregate sayılar düzeltildi.

## Henüz yapılmadı / sıradaki somut adımlar

1. **Üç modelin karşılaştırılması**: B. subtilis (su-kısıtlı uçurum, Mars
   esansiyelliği DEĞİŞTİRMİYOR), Salinibacter (uçurum yok, doğrusal,
   esansiyellik değişmiyor), JCVI-syn3A (glikoz-kısıtlı uçurum, Mars 4
   YENİ esansiyel gen EKLİYOR — PDH/PTA/ACK yolu). Bu üç farklı
   "kısıtlayıcı darboğaz + esansiyellik tepkisi" profili muhtemelen
   makalenin ana karşılaştırma bulgusu olabilir; ayrı bir
   `karsilastirma.py`/tablo ile bir araya getirilmeli.
2. **README'de dürüstlük notu eklenmeli**: referans ortamın (zengin/
   tanımsız) B. subtilis/Salinibacter'deki tanımlı-ortam referanslarından
   niteliksel farkı, makalede/README'de MUTLAKA vurgulanmalı — aksi halde
   üç model arası büyüme oranı karşılaştırması yanıltıcı olur.
3. ATPase-bağlı 8 genin essentiality paradoksu için literatürde JCVI-
   syn3A'da oksidatif fosforilasyonun/ATP sentazın gerçek rolü araştırılıp
   README'ye bir yorum eklenebilir (şu an sadece model sınırlaması olarak
   not düşüldü, biyolojik yorum yapılmadı).

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
