"""
Mars kısıt şiddetine duyarlılık analizi (sensitivity analysis) — JCVI-syn3A.

mars-minimal-gene-network/src/mars_duyarlilik.py ile AYNI yöntem ve AYNI
şiddet-ekseni ankorları (SERT/ILIMLI) kullanılıyor -- bu, iki proje
arasında doğrudan karşılaştırılabilir bir t-ekseni sağlıyor. Nedenler
(radyasyonun ATPM'ye sayısal etkisinin bilinmemesi, atmosfer yüzdelerini
akı sınırına çevirecek kinetik veri eksikliği) için bkz. o dosyanın
docstring'i / bu projenin README'si.

İlk (README'deki) ön bulgu: burada baskın kısıt B. subtilis'teki gibi su
değil, GLİKOZ -- keskin bir uçurum glc_lb ≈ -0.8/-0.75 mmol/gDW/h
arasında bulundu (o2/h2o tek başına belirleyici değil). Bu script, aynı
üç-parametreli birleşik tarama ile bu bulguyu t-ekseni üzerinde
doğruluyor/görselleştiriyor.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from mars_fba import bakim_reaksiyonunu_bul, mars_kisitlarini_uygula, modeli_yukle, referans_buyume

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SONUC_KLASORU = os.path.join(PROJE_KOKU, "results")

# mars-minimal-gene-network ile AYNI ankorlar (bkz. o projenin docstring'i).
SERT = dict(o2=-0.5, glc=-0.05, h2o=1.0)
ILIMLI = dict(o2=-10.0, glc=-2.0, h2o=30.0)

BAKIM_CARPANLARI = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
T_DEGERLERI = [round(i * 0.02, 2) for i in range(51)]  # 0.00, 0.02, ..., 1.00


def kisit_degerleri(t):
    o2 = SERT["o2"] + t * (ILIMLI["o2"] - SERT["o2"])
    glc = SERT["glc"] + t * (ILIMLI["glc"] - SERT["glc"])
    h2o = SERT["h2o"] + t * (ILIMLI["h2o"] - SERT["h2o"])
    return o2, glc, h2o


def tek_nokta_calistir(model, atpase, atpase_taban, t, bakim_carpani):
    o2, glc, h2o = kisit_degerleri(t)
    # bakim_taban'ı MUTLAKA açıkça veriyoruz (bkz. mars_fba.py uyarısı):
    # model/atpase bu döngü boyunca tekrar kullanılıyor, atpase.upper_bound
    # bir önceki çağrıdan kalma mutasyona uğramış bir değer olabilir.
    mars_kisitlarini_uygula(
        model, atpase, o2_lb=o2, glc_lb=glc, h2o_cap=h2o, bakim_carpani=bakim_carpani,
        bakim_taban=atpase_taban, sessiz=True,
    )
    sol = model.optimize(raise_error=False)
    durum = model.solver.status
    buyume = sol.objective_value if durum == "optimal" else None
    return durum, buyume, o2, glc, h2o


def tarama_yap(model, atpase, atpase_taban):
    satirlar = []
    for bakim_x in BAKIM_CARPANLARI:
        for t in T_DEGERLERI:
            durum, buyume, o2, glc, h2o = tek_nokta_calistir(model, atpase, atpase_taban, t, bakim_x)
            satirlar.append(
                dict(
                    bakim_carpani=bakim_x,
                    t=t,
                    O2_lb=o2,
                    glc_lb=glc,
                    h2o_cap=h2o,
                    durum=durum,
                    buyume=buyume,
                )
            )
    return pd.DataFrame(satirlar)


def grafik_ciz(df, baseline_buyume, dosya_yolu):
    fig, ax = plt.subplots(figsize=(8, 5))
    for bakim_x, grup in df.groupby("bakim_carpani"):
        gecerli = grup.dropna(subset=["buyume"]).sort_values("t")
        ax.plot(gecerli["t"], 100 * gecerli["buyume"] / baseline_buyume, marker=".", label=f"bakım ×{bakim_x}")
    ax.set_xlabel("Şiddet ekseni t  (0 = en sert ilk varsayım, 1 = ılımlı uç)")
    ax.set_ylabel("Büyüme oranı (referansa göre %)")
    ax.set_title("Mars kısıt şiddetine duyarlılık analizi — JCVI-syn3A (iMMSYN)")
    ax.legend(title="Bakım çarpanı")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(dosya_yolu, dpi=150)
    plt.close(fig)
    print(f"Grafik kaydedildi: {dosya_yolu}")


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    model = modeli_yukle()
    baseline = referans_buyume(model)  # kısıtlar uygulanmadan ÖNCE ölçülmeli
    atpase = bakim_reaksiyonunu_bul(model)
    atpase_taban = atpase.lower_bound  # standart yöne çevrilmiş orijinal 0.575 -- döngü boyunca hep bundan çarpıyoruz

    df = tarama_yap(model, atpase, atpase_taban)

    csv_yolu = os.path.join(SONUC_KLASORU, "duyarlilik_sonuclari.csv")
    df.to_csv(csv_yolu, index=False)
    print(f"\nSonuçlar kaydedildi: {csv_yolu} ({len(df)} satır)")

    grafik_yolu = os.path.join(SONUC_KLASORU, "buyume_vs_siddet.png")
    grafik_ciz(df, baseline.objective_value, grafik_yolu)

    print("\n--- Özet: her bakım çarpanı için ilk feasible nokta ---")
    for bakim_x, grup in df.groupby("bakim_carpani"):
        gecerli = grup.dropna(subset=["buyume"]).sort_values("t")
        if gecerli.empty:
            print(f"bakım ×{bakim_x}: taranan aralıkta hiçbir t değeri feasible değil")
            continue
        ilk = gecerli.iloc[0]
        pct = 100 * ilk["buyume"] / baseline.objective_value
        print(
            f"bakım ×{bakim_x}: ilk feasible t={ilk['t']:.2f} "
            f"(O2={ilk['O2_lb']:.2f}, glc={ilk['glc_lb']:.2f}, h2o=±{ilk['h2o_cap']:.2f}) "
            f"-> büyüme={ilk['buyume']:.6f} (%{pct:.1f} referans)"
        )


if __name__ == "__main__":
    main()
