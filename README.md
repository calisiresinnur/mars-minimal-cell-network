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

**Henüz yapılmadı**: tam duyarlılık analizi (B. subtilis/Salinibacter
projelerindeki `mars_duyarlilik.py` gibi, tüm parametre kombinasyonlarını
tarayan), gen esansiyellik/silme analizi (`SOLVER_TOLERANCE=1e-9` ile),
ve bu bulguların B. subtilis/Salinibacter karşılaştırmasıyla bir araya
getirilmesi.

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
    └── mars_fba.py           # Model yükleme + Mars kısıtları + FBA
```

## Kaynaklar

- Breuer M, Käser T, Kolar K, et al. (2019) Essential metabolism for a
  minimal cell. *eLife* 8:e36842. <https://doi.org/10.7554/eLife.36842>
  (PMID: 30657448, PMC6609329)
- Diğer atmosfer/radyasyon/NGAM kaynakları için bkz.
  [mars-minimal-gene-network README](https://github.com/calisiresinnur/mars-minimal-gene-network) —
  aynı Mars koşul kısıtları (O2/CO2/su/radyasyon) burada da geçerli.
