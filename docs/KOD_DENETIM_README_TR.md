# Alzheimer three-omics source-code audit — Türkçe

Bu klasördeki Archived_Original_Analysis_Modules, önceki çalışma arşivlerinden **değiştirilmeden** alınmış özgün analiz betikleridir. Her modül tek başına tüm projeyi üretmez; varsayılan dosya yolları, internet erişimi, gerekli R/Python paketleri ve özgün girdi dosyaları ayrı ayrı sağlanmalıdır.

## Analiz zinciri (orijinal modül adları)

06: GSE33000 iki renkli microarray hazırlama ve yaş/cinsiyet ayarlı `limma` modeli (R+limma). 07: GPL4372 prob–gen eşleme ve gen düzeyi özetleme. 10: PXD010560 yayımlanmış destek tablolarını alma, protein/transkript eşleştirme. 17: hedef CpG kovaryat duyarlılığı. 18: GSE134379 MTG tam CpG düzeyi yaş/cinsiyet/plate ayarlı OLS. 23: anotasyon tabanlı SNP-risk dışlama ve BH-FDR tekrar hesabı. 24: SNP-risk filtreli `missMethyl`. 25: üç-omik betimleyici entegrasyon.

Mevcut `verify_supplementary.py` yalnızca ek tablolardan sayı tutarlılığını yeniden denetler; `run_gometh.R` yalnızca hazır CpG listeleriyle enrichment hesaplar. `python_validation_retest.txt` ve `code_audit.json` gerçek test sonuçlarını içerir.

## Eksikler ve yeniden üretilebilirlik sınırı

- **Ham veriden bütün zincir bu ortamda çalıştırılmadı.** Özellikle GSE134379 MTG AD ham beta matrisi sadece `.part` dosyası olarak bulundu; tam arşiv yok. Bu dosya tamamlanmadan Module18 sonucu bağımsız yeniden hesaplanamaz.
- PXD010560 kaynak destek tabloları, GSE33000 GEO ham dosyaları ve gerekli ara çıktıların modüllerin beklediği dizin yapısına uygun kurulması gerekir.
- R çalışma ortamı ve Bioconductor bağımlılıkları bulunmuyorsa `limma` ve `missMethyl` gerçek çalışma testi mümkün değildir.
- İçerilen tüm Python betikleri sözdizimi açısından kontrol edilmiştir; bu **gerçek veriyle uçtan uca doğrulama değildir**.
- `run_gometh.R` bazı modüllerle aynı işleve yöneliktir; birden fazla script tek bir doğrulanmış sürüm olarak birleştirilmemiştir.
- Ek tabloların kontrolü, ham veri kalitesini veya bilimsel sonuçların bağımsız doğruluğunu kanıtlamaz.

## Çalıştırma

Önce `python Code/verify_supplementary.py --folder . --write-gometh-inputs Code/Module21_inputs`.
Sonra uygun R ortamında: `Rscript Code/run_gometh.R --inputs Code/Module21_inputs --out Code/Module21_results --directions`.
Ham veri analizi için önce her modülün girdi gereksinimlerini ve sürümlerini kontrol edin; eksik ham dosyalara rağmen final sayıları yeniden üretilmiş gibi sunmayın.
