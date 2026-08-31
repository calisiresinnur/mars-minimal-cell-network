# mars-minimal-cell-network

Mars yüzey koşullarında bir **minimal sentetik hücrenin** (JCVI-syn3A,
genom-ölçekli metabolik model + FBA ile) metabolik olarak canlı kalıp
kalamayacağını hesaplamalı incelemek.

Bu proje, [mars-minimal-gene-network](https://github.com/calisiresinnur/mars-minimal-gene-network)
projesinin devamı/kardeş projesidir. O projede "Dünya'da yaşayan doğal bir
bakterinin (B. subtilis) minimal gen seti Mars'ta da yeterli mi" sorusu
soruldu; burada soru tersine çevriliyor: **"Mars'ta hayatta kalmak için
teorik olarak en az kaç/hangi gen yeterli?"** — bunun için doğal bir
organizma yerine, zaten deneysel olarak minimuma indirgenmiş sentetik bir
genomu (JCVI-syn3A) başlangıç noktası alıyoruz.

IAC 2026 · IAF/IAA Space Life Sciences Symposium (A1), Paper ID 114761
kapsamındaki araştırmanın bir parçasıdır. Yazar: Esinnur Çalışır,
İstanbul Üniversitesi.

## Model

**iMMSYN** — Breuer, Käser, Kolar, Alcaraz, Grote, Strauss, Vashee, Suthers,
Vitkin, Peckham, Fraser, Kruse, Smith, Glass, Palsson, Elhai, Baliga, 2019,
eLife, *"Essential metabolism for a minimal cell"*, DOI:
[10.7554/eLife.36842](https://doi.org/10.7554/eLife.36842) (PMID: 30657448,
PMC6609329).

JCVI-syn3A için: **155 gen, 338 reaksiyon, 304 metabolit**. Deneysel
transpozon mutajenez verisiyle doğrulanmış (in vivo esansiyellik %92, in
silico %79, Matthews correlation 0.59).

Model dosyası makalenin **Supplementary file 9**'u (SBML/FBC formatı),
eLife'ın kendi "Figures and data" sayfasından doğrulanarak indirildi:
<https://cdn.elifesciences.org/articles/36842/elife-36842-supp9-v3.xml.zip>
→ yerelde `data/models/iMMSYN.xml.gz` olarak önbelleğe alındı.

**Doğrulama**: `cobra.io.read_sbml_model` ile yüklendiğinde model
355/338/304 gen/reaksiyon/metabolit sayısını makaleyle birebir veriyor;
yayınlanmış (zengin/tanımsız) varsayılan ortamda büyüme oranı **0.342/saat**
hesaplandı — makalenin bildirdiği ~2 saatlik çiftlenme süresiyle
(ln2/2 ≈ 0.347/saat) örtüşüyor.

### mars-minimal-gene-network'e göre önemli yapısal farklar

- **Referans ortam**: iYO844 (B. subtilis) ve iMB631 (Salinibacter) için
  tanımlı/kalibre edilmiş minimal ortamlar kullanılmıştı. JCVI-syn3A için
  literatürde "normal büyümeyi destekleyen tanımlı bir ortam henüz elde
  edilmedi" (Breuer ve ark. 2019) — model, yayınlanmış haliyle
  zengin/tanımsız bir ortam varsayıyor (78 exchange bileşeninin çoğu
  sınırsız alım, glikoz tek enerji kaynağı). Bu proje, makalenin kendi
  metodolojisiyle tutarlı olarak bu varsayılan ortamı "Dünya benzeri"
  referans olarak kullanıyor — bu, önceki iki modeldeki tanımlı-ortam
  referanslarından **niteliksel olarak farklı bir varsayım**.
- **Bakım enerjisi (NGAM)**: iYO844'teki gibi ayrı, sabit (lb=ub) bir
  "ATPM" reaksiyonu YOK. Makale (Bölüm 'GAM/NGAM') NGAM'ı üç bileşene
  dağıtıyor: (1) `ATPase` reaksiyonu üzerinde tek yönlü bir sınır (makalede
  0.57 mmol/gDW/h — modelde tam olarak **0.575**, birebir doğrulandı),
  (2) `Protein_degrad` alt sınırı (makalede 3.5×10⁻⁴ — modelde **0.00035**,
  birebir), (3) `RNA_degrad` alt sınırı (makalede 7.7×10⁻³ — modelde
  **0.0077**, birebir). Bu proje, önceki projedeki "bakım çarpanı"
  yöntemini `ATPase` reaksiyonunun sabit akışına uygulayarak genişletiyor.
  NOT: `ATPase` reaksiyonu modelde ters yönde yazılmış olduğu için
  sabitleyen kısıt `lower_bound` değil `upper_bound`'dur (bkz.
  `src/mars_fba.py` docstring'i).

## Yöntem: Windows Unicode kullanıcı adı sorunu

`mars-minimal-gene-network`'te keşfedilen aynı çözüm burada da geçerli:
Windows'ta kullanıcı adında ASCII-dışı karakter ("Ergün") olduğunda
libSBML dosya yolunu C seviyesinde açamıyor. Bunu aşmak için `.xml.gz`
dosyası Python'ın gzip modülüyle açılıp ham SBML metni libSBML'e dosya
yolu değil **string** olarak veriliyor.

## Şu ana kadarki bulgu (ilk geçiş, henüz duyarlılık analizi yapılmadı)

`src/mars_fba.py`, önceki projedeki B. subtilis parametreleriyle
birebir aynı başlangıç değerleriyle (O2 alımı ≥ -0.5, glikoz alımı ≥ -0.05,
su akışı ±1.0, bakım çarpanı ×3) çalıştırıldığında model **infeasible**
çıkıyor.

Kısıtları teker teker izole ederek (diğerlerini tamamen açık bırakarak)
yapılan bir ön tarama, B. subtilis'teki bulgudan **niteliksel olarak
farklı** bir sonuç veriyor:

- B. subtilis (iYO844)'te asıl kısıtlayıcı su akışıydı (`h2o_cap`).
- JCVI-syn3A'da su akışı tek başına kısıtlayıcı DEĞİL (su kısıtı=1.0
  tek başına feasible, büyüme=0.328/saat, referansa çok yakın).
- Bunun yerine **glikoz alımı** baskın kısıt: `EX_glc__D_e` alt sınırı
  ile keskin bir feasibility uçurumu var, eşik **-0.8 ile -0.75
  mmol/gDW/h arasında** (diğer kısıtlar tamamen açıkken): -0.8'de hâlâ
  feasible (büyüme ≈0.0016/saat), -0.75'te tamamen infeasible.
- Bu, biyolojik olarak beklenir: JCVI-syn3A zaten genomu en aza
  indirilmiş, alternatif karbon/enerji kaynağı yollarından büyük
  ölçüde arındırılmış bir organizma — glikoz kesildiğinde yedek bir
  metabolik rota yok.

**Önemli metodolojik not (bu oturumda canlı yakalanan bir solver
artefaktı)**: `model.optimize()` art arda birden çok kez (önce referans,
sonra Mars kısıtlarıyla) çağrıldığında, glpk/optlang infeasible bir
çözüm için bazen sıfır olmayan, anlamsız bir `objective_value`
döndürebiliyor (warm-start solver durumunun taşınması nedeniyle).
`src/mars_fba.py`, `model.solver.status`'u EXPLICIT olarak kontrol
ediyor ve infeasible durumda `objective_value`'yu asla gerçek bir
büyüme oranı gibi raporlamıyor. Bu, `mars-minimal-gene-network`'teki
tolerance-artefaktı dersiyle aynı ailede bir uyarı — SOLVER_TOLERANCE
zaten 1e-9'a çekili (bkz. `src/mars_fba.py`).

## Duyarlılık analizi

`src/mars_duyarlilik.py`, B. subtilis/Salinibacter projeleriyle AYNI
şiddet-ekseni (t: 0=sert, 1=ılımlı) ve aynı bakım-çarpanı listesini
kullanıyor. Sonuç, yukarıdaki ön bulguyu doğruluyor: keskin feasibility
uçurumu, bakım çarpanına göre t≈0.40 (çarpan×1.0) ile t≈0.62 (çarpan×4.0)
arasında. **Bu modelde uçurum, B. subtilis'tekinden bile daha keskin** —
t*'yi 6 ondalık basamağa yuvarlayıp elle kopyalamak bile (~1e-6 hassasiyet
kaybı) "optimal"ı "infeasible"a çevirebiliyor (bkz. `results/`).
Sonuçlar: `results/duyarlilik_sonuclari.csv`, `results/buyume_vs_siddet.png`.

## Gen esansiyellik/silme analizi — DOĞRULANMIŞ SONUÇ

`src/mars_gen_silme.py`, referans (kısıtsız) + üç Mars senaryosu
(bakım çarpanı ×1.5/×2.0/×3.0, her biri WT büyümenin referansın en az
%5'i olduğu, 5x bağımsız tekrarla doğrulanmış "rahat" bir noktada) için
155 genin tek tek silinmesini test ediyor.

**Bu analiz sırasında ÜÇ ayrı hata canlı yakalanıp düzeltildi** (ayrıntı
için `src/mars_gen_silme.py` docstring'i) — kullanıcının isteğiyle
yapılan bir bilimsel doğruluk denetimi sırasında, sadece sayısal değil
mantıksal olarak da:

1. cobra'nın standart gen-silme mekanizması, infeasible olan KO'larda
   `growth` alanını `NaN` bırakıyor; `oran = growth/wt_büyüme` de `NaN`
   olunca `oran < eşik` pandas'ta HER ZAMAN `False` dönüyor — yani
   infeasible (en kesin esansiyel durum!) YANLIŞLIKLA "esansiyel değil"
   sayılıyordu. İlk çalıştırmada bu, "0 esansiyel gen" gibi biyolojik
   olarak imkânsız bir sonuca yol açtı.
2. NGAM'ın `ATPase` bileşenine bağlı 8 gen (F1F0-ATP sentaz alt
   birimleri, MMSYN1_0789-0796) VE `Protein_degrad`'a bağlı 1 gen
   (MMSYN1_0394) — toplam 9 gen — "silindiğinde" büyüme ARTIYOR
   (bounds (0,0)'a çekiliyor, bu kapasiteyi değil zorunlu NGAM YÜKÜNÜ
   kaldırıyor). İlk yorumum "model sınırlaması, hata değil" şeklindeydi
   — ama makalenin kendi **Table 4**'ünü (locus-düzeyinde Ess_FBA
   sütunu) XML tam metninden çekip çapraz kontrol ettiğimde, makalenin
   BU 9 GENİN TAMAMINI esansiyel (■) işaretlediğini gördüm. Yani sorun
   modelin yapısı değil, benim naif "bounds->(0,0)" knockout
   yorumumdu: NGAM'ın dayattığı minimum akış hücrenin gerçek, sürekli
   bir fizyolojik ihtiyacını temsil ediyor — ilgili gen silindiğinde bu
   ihtiyaç ORTADAN KALKMAZ, sadece KARŞILANAMAZ hale gelir (infeasible
   olmalı). Düzeltme: bu 9 gen artık ham simülasyon sonucundan bağımsız
   olarak esansiyel kabul ediliyor.
3. Bu iki düzeltme sonrası **referans esansiyel gen sayısı 123/155
   (%79.4)** çıktı — makalenin kendi bildirdiği 123/155 (%79) ile
   örtüşüyor (bkz. makale metni: *"123 of the 155 genes included in the
   model are essential (79%)"*). Ayrıca birkaç bağımsız gen (5 tRNA
   sentetaz: ALATRS/0163, ARGTRS/0535, ASNTRS/0076, ASPTRS/0287,
   CYSTRS/0837 — makalede Ess_FBA=■; ve 3 non-essential kontrol geni:
   MMSYN1_0330/0382/0876 — makalede işaretsiz/non-essential) tek tek
   makaleyle karşılaştırıldı, hepsi örtüştü.

   **Metodolojik dürüstlük notu (kullanıcının uyarısıyla eklendi)**: madde
   2'deki düzeltmeye giden yol şöyle işledi — önce 114+9=123 sayısal
   eşleşmesini fark ettim, SONRA makale Table 4'e bakıp bu 9 genin orada
   esansiyel işaretli olduğunu gördüm. Yani "önce sayı, sonra gerekçe"
   sırası, gerekçeyi sayıya uydurma riski taşıyor — bunu saklamıyorum.
   Gerekçenin (NGAM = gerçek, süregelen fizyolojik ihtiyaç; gen silindiğinde
   ihtiyaç değil karşılanabilirlik ortadan kalkar; standart GEM
   essentiality konvansiyonu bu tür genleri esansiyel sayar) sayıdan
   BAĞIMSIZ olarak da savunulabilir olduğuna inanıyorum, ama bu benim
   değerlendirmem — okuyucu ham veriyi (`esansiyel_ham` sütunu,
   düzeltmesiz: referansta 114/155) görüp kendi yargısını verebilir. Bu
   satır, referans senaryo için geçerli (makalede karşılaştıracak veri
   olduğu için). **Mars'a özgü asıl bulgu (aşağıdaki 4 gen) bu düzeltmeden
   TAMAMEN BAĞIMSIZ** — makalede Mars kısıtları altında essentiality
   diye bir veri yok, dolayısıyla o bulguda "makaleye uydurma" riski
   yapısal olarak söz konusu değil; düzeltme öncesi de sonrası da aynı
   4 gen çıkıyor.
4. Ayrıca makalenin kendi metni, referans (Dünya benzeri) koşulda tam
   olarak benim bulduğum 4 genin (aşağıya bkz.) doğrulanma süresini
   (doubling time) ayrıca not düşmüş: *"single knockouts of loci
   pdhC/0227 through ackA/0230 ... had doubling times of 3.22 hr"*
   (referans 2.02 hr'a karşı). Oran: 2.02/3.22=0.6273 — benim
   hesapladığım oran (0.6289) ile ~%0.2 farkla örtüşüyor. Bağımsız bir
   ikinci doğrulama.

**Ana bulgu — Mars'a özgü 4 yeni esansiyel gen** (üç Mars senaryosunun
DA HEPSİNDE tutarlı, referansta esansiyel değil, yukarıdaki düzeltmelerden
etkilenmiyor): `MMSYN1_0227 (pdhC), MMSYN1_0228, MMSYN1_0229 (pta),
MMSYN1_0230 (ackA)`. Bu dört gen birlikte **piruvat dehidrogenaz (PDH) →
fosfotransasetilaz (PTAr) → asetat kinaz (ACKr)** yolunu kodluyor —
piruvattan asetil-CoA/NADH üretimi, ardından asetat kinaz üzerinden
**ek ATP** (substrat düzeyinde fosforilasyon) üretimi. Yorum: zengin/
kısıtsız ortamda hücre bu ekstra ATP kaynağı olmadan da (daha yavaş,
makalenin kendi doğruladığı gibi) büyüyebiliyor; ama Mars'ın
sıkılaştırılmış enerji bütçesinde (kısıtlı glikoz/O2 + artmış bakım) bu
yol tamamen vazgeçilmez hale geliyor (referansta ölümcül değil, Mars'ta
ölümcül). Bu, B. subtilis projesindeki `dltABCD` bulgusunun (Mars'ta
esansiyellikten ÇIKAN genler) TERSİ yönünde ama aynı ailede bir bulgu:
enerji darboğazı, ATP üretimine katkısı olan HER yolu kritikleştiriyor.

**Sonuç tablosu (düzeltilmiş)**:

| Senaryo | Esansiyel gen | Oran |
|---|---|---|
| Referans (kısıtsız) | 123/155 | %79.4 (makale: %79) |
| Mars ×1.5 / ×2.0 / ×3.0 (üçü de) | 127/155 | %81.9 |

Sonuçlar: `results/gen_silme_sonuclari.csv` (hem ham `esansiyel_ham` hem
düzeltilmiş `esansiyel` sütunları içerir, şeffaflık için),
`results/mars_yeni_esansiyel_genler.csv`,
`results/mars_dispanse_olan_genler.csv` (bu proje için boş — hiçbir gen
dispanse olmuyor).

**Henüz yapılmadı**: B. subtilis/Salinibacter/JCVI-syn3A karşılaştırmalı
"kısıtlayıcı darboğaz" tablosunun bir araya getirilmesi.

## Kurulum

```bash
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
python src/mars_fba.py
```

## Repo yapısı

```
.
├── README.md
├── DEVAM_NOTLARI.md
├── requirements.txt
├── data/
│   └── models/
│       └── iMMSYN.xml.gz     # JCVI-syn3A modeli (Breuer ve ark. 2019, Supp. file 9)
├── results/
└── src/
    ├── mars_fba.py           # Model yükleme + Mars kısıtları + FBA
    ├── mars_duyarlilik.py    # Duyarlılık analizi + grafik
    └── mars_gen_silme.py     # Tekli gen silme (SOLVER_TOLERANCE=1e-9)
```

## Kaynaklar

- Breuer M, Käser T, Kolar K, et al. (2019) Essential metabolism for a
  minimal cell. *eLife* 8:e36842. <https://doi.org/10.7554/eLife.36842>
  (PMID: 30657448, PMC6609329)
- Diğer atmosfer/radyasyon/NGAM kaynakları için bkz.
  [mars-minimal-gene-network README](https://github.com/calisiresinnur/mars-minimal-gene-network) —
  aynı Mars koşul kısıtları (O2/CO2/su/radyasyon) burada da geçerli.
