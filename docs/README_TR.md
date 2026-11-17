# Alzheimer çoklu-omik çalışması — yeniden üretilebilirlik dosyaları

## Dosyalar
- `verify_supplementary.py`: Python 3 standart kütüphanesiyle 11.471 aday CpG, 7.382/4.089 metilasyon yönü, altı GO/KEGG zenginleşme sayısı, 74 protein–transkript eşleştirmesi ve 60/14 yön uyumunu orijinal supplementary CSV dosyaları üzerinden yeniden hesaplar. Ayrıca R enrichment betiği için dört giriş listesi üretir.
- `run_gometh.R`: Önceki Module21 çalışmasından korunan missMethyl `gometh` analizi; `all.cpg` arka planı, `prior.prob=TRUE` ve GO/KEGG analizini kullanır. **Bu betik ham methylation regresyonunu yeniden üretmez**; SNP-risk dışlama işlemi betik içinde değil, sunulan filtrelenmiş CpG girdi listelerinde uygulanır.
- `Module21_inputs/`: Sağlanan supplementary S5 ve S6 CSV dosyalarından deterministik olarak yeniden üretilmiş filtrelenmiş tested/candidate/hypo/hyper CpG listeleri.
- `validation_output.txt`: Python betiğinin test çıktısı.

## Kullanım
ZIP paketi aynı dizine açtıktan sonra:

`python Code/verify_supplementary.py --folder . --write-gometh-inputs Code/Module21_inputs`

R/Bioconductor'da `missMethyl`, `jsonlite` ve `IlluminaHumanMethylation450kanno.ilmn12.hg19` kuruluysa:

`Rscript Code/run_gometh.R --inputs Code/Module21_inputs --out Code/Module21_results --directions`

## Şeffaflık
Python betiği **özet sonuçların ek tablolarla tutarlılığını** kontrol eder; GEO/PRIDE ham veri indirme, QC, GSE33000 limma modeli, GSE134379 yaş/cinsiyet/plate ayarlı OLS ve proteomik özgün karşılaştırmasını uçtan uca **yeniden üretmez**. Orijinal bu basamakların tekil yürütülebilir betikleri bu pakette yoktur. R analizinin gerçek çalıştırılması bu ortamda doğrulanmamıştır. Bu kapsam sınırlaması makalenin Code availability bölümünde açıklanmıştır.
