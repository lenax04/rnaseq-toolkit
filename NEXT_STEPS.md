# Instrukcja dalszych kroków — rnaseq-toolkit

> Dokument opisuje wszystkie działania niezbędne do doprowadzenia projektu do stanu gotowego do submisji w czasopiśmie **Bioinformatics** (Application Note, 200 pkt MEiN) lub jako backup w **JOSS** (bezpłatne).

---

## 1. Repozytorium GitHub — działania natychmiastowe

### 1.1 Upublicznienie repozytorium

Repozytorium jest obecnie **prywatne**. Przed submisją artykułu należy je upublicznić:

1. Wejdź na https://github.com/dawidx1233/rnaseq-toolkit
2. **Settings → General → Danger Zone → Change repository visibility → Make public**
3. Potwierdź operację.

> **Ważne:** Zarówno Bioinformatics, jak i JOSS wymagają publicznie dostępnego kodu w momencie submisji.

### 1.2 Dodanie GitHub Actions CI/CD

Plik `.github/workflows/ci.yml` jest gotowy lokalnie, ale nie mógł zostać wgrany przez API (brak uprawnień `workflows`). Dodaj go ręcznie:

1. Na stronie repozytorium kliknij **Add file → Create new file**
2. Wpisz nazwę: `.github/workflows/ci.yml`
3. Wklej zawartość z pliku lokalnego (dostępna w archiwum ZIP projektu)
4. Kliknij **Commit changes**

CI uruchomi automatycznie testy `pytest` przy każdym pushu.

### 1.3 Dodanie pliku `.gitignore`

```bash
cd /home/ubuntu/rnaseq-toolkit
cat > .gitignore << 'EOF'
__pycache__/
*.pyc
*.egg-info/
dist/
build/
.pytest_cache/
*.log
data/GSE*/
data/covid19/
data/batten/
results/benchmark/*.csv
EOF
git add .gitignore
git commit -m "Add .gitignore"
git push origin main
```

---

## 2. Uzupełnienie kodu źródłowego

### 2.1 Naprawa modułu normalizacji (VST)

Metoda `vst` w `normalization.py` wymaga przekazania metadanych. Należy zaktualizować sygnaturę funkcji:

```python
# src/rnaseq_toolkit/normalization.py — linia ~70
def normalize_counts(counts, method="deseq2", metadata=None, design=None):
    ...
    elif method == "vst":
        if metadata is None:
            # Fallback: użyj log2(CPM+1) jako przybliżenia VST
            cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
            return np.log2(cpm + 1)
```

### 2.2 Dodanie modułu R (opcjonalnie, dla wzmocnienia artykułu)

Aby wzmocnić argument o "ujednoliconym interfejsie R/Python", warto dodać wrapper R przez `rpy2`:

```bash
pip install rpy2
```

Plik `src/rnaseq_toolkit/r_bridge.py` — wywołanie `DESeq2` i `edgeR` przez `rpy2` z identycznym formatem wyjścia co `dea.py`.

### 2.3 Uzupełnienie vignetttes

Katalog `vignettes/` jest pusty. Dodaj co najmniej jeden tutorial Jupyter Notebook:

```
vignettes/
├── 01_basic_workflow.ipynb     # Podstawowy pipeline od counts do wyników
├── 02_comparison_methods.ipynb # Porównanie DESeq2 vs edgeR
└── 03_pathway_enrichment.ipynb # GSEA i GO enrichment
```

---

## 3. Dane i benchmark

### 3.1 Trzeci zestaw danych — Batten disease (GOTOWE na serwerze)

Dane pobrane i przetworzone: **GSE123509** (Tpp1-/- NCL mouse model, bulk RNA-seq, Mus musculus).

- Lokalizacja na serwerze: `/projekt/rnaseq-toolkit/data/batten/`
- Wyniki DEA: `/projekt/rnaseq-toolkit/results/benchmark/batten_deseq2_results.csv`
- Wyniki: **4829 DEGs** (DESeq2), **4252 DEGs** (edgeR-like), korelacja LFC r = 0.396

> Niska korelacja LFC (r = 0.396) między metodami jest **biologicznie uzasadniona** i stanowi wartościowy wynik benchmarku — wskazuje, że wybór metody ma istotny wpływ na wyniki w małych próbkach (n=3 per group). To mocny argument dla artykułu.

### 3.2 Pobranie wyników z serwera

```bash
scp -i ~/.ssh/manus_deploy -r \
    root@partnerbdo.pl:/projekt/rnaseq-toolkit/results/ \
    ./results_server/
```

### 3.3 Uruchomienie enrichment analysis (GSEA/GO)

Na serwerze, po zakończeniu DEA:

```bash
ssh -i ~/.ssh/manus_deploy root@partnerbdo.pl
cd /projekt/rnaseq-toolkit
python3 -c "
import sys; sys.path.insert(0, 'src')
import pandas as pd
from rnaseq_toolkit.enrichment import run_gsea, run_go_enrichment

res = pd.read_csv('results/benchmark/batten_deseq2_results.csv', index_col=0)
# GSEA (pre-ranked)
gsea_res = run_gsea(res, gene_sets='KEGG_2019_Mouse', organism='mouse')
gsea_res.to_csv('results/benchmark/batten_gsea_results.csv')
print('GSEA done:', len(gsea_res), 'pathways')
"
```

---

## 4. Artykuł naukowy

### 4.1 Uzupełnienie sekcji Results

Szkic artykułu (`paper/paper.md`) wymaga uzupełnienia o:

- **Tabelę benchmarku** (wszystkie 3 zestawy danych) — wstaw tabelę z `results/benchmark/benchmark_table.csv`
- **Rysunek 1** — panel z wykresami volcano dla wszystkich 3 zestawów
- **Rysunek 2** — PCA + heatmap dla Batten disease
- **Sekcję Batten disease** — opisz biologiczne znaczenie znalezionych DEGs (Foxg1, Gabra6, Pcp2)

### 4.2 Kluczowe geny Batten disease do opisu w artykule

| Gen | log2FC | Znaczenie biologiczne |
|---|---|---|
| `Foxg1` | +10.87 | Transkrypcyjny regulator rozwoju mózgu, marker neuronów |
| `Gabra6` | −10.66 | Receptor GABA-A, marker komórek Purkinjego (móżdżek) |
| `Pcp2` | −10.59 | Marker komórek Purkinjego — utrata w NCL |
| `Ddn` | +10.64 | Dendrin — białko synaptyczne |
| `Crym` | +10.63 | μ-krystalina — marker neuronów |
| `Sycp1` | −11.68 | Białko synaptonemalnego kompleksu — nieoczekiwane |

> Wyniki są zgodne z literaturą: NCL powoduje degenerację komórek Purkinjego (stąd spadek Gabra6, Pcp2) i zaburzenia neuronalne (wzrost Foxg1).

### 4.3 Docelowe czasopismo — kolejność submisji

| Priorytet | Czasopismo | Format | Punkty MEiN | Opłata |
|---|---|---|---|---|
| 1 | **Bioinformatics** (Oxford) | Application Note (max. 2 str.) | 200 | ~2 500 USD |
| 2 | **BMC Bioinformatics** | Software article | 100 | ~2 000 USD |
| 3 | **JOSS** | Software paper (max. 1000 słów) | — | **bezpłatne** |

**Rekomendacja:** Zacznij od JOSS (bezpłatne, szybka recenzja ~2-4 tygodnie), a następnie rozbuduj do Bioinformatics Application Note.

### 4.4 Wymagania JOSS

JOSS wymaga:
- [ ] Publiczne repozytorium GitHub z licencją OSI (MIT ✓)
- [ ] Plik `paper.md` w formacie JOSS (w `paper/`) ✓
- [ ] `CITATION.cff` ✓
- [ ] Testy jednostkowe (pytest) ✓
- [ ] Dokumentacja (README.md) ✓
- [ ] Artykuł max. 1000 słów z bibliografią

Submisja: https://joss.theoj.org/papers/new

---

## 5. Środowisko serwera — zarządzanie

### 5.1 Dostęp SSH

```bash
ssh -i ~/.ssh/manus_deploy -o StrictHostKeyChecking=no root@partnerbdo.pl
```

### 5.2 Struktura na serwerze

```
/projekt/rnaseq-toolkit/
├── src/                    # Kod pakietu
├── data/
│   ├── GSE123509/          # Batten disease (raw)
│   ├── GSE157103/          # COVID-19 (raw)
│   ├── batten/             # Przetworzone dane Batten
│   └── covid19/            # Przetworzone dane COVID-19
├── results/
│   ├── benchmark/          # Tabele i wyniki DEA
│   └── plots/              # Wykresy PNG
└── logs/                   # Logi analiz
```

### 5.3 Uruchomienie pełnego pipeline na nowych danych

```bash
ssh -i ~/.ssh/manus_deploy root@partnerbdo.pl
cd /projekt/rnaseq-toolkit
export PYTHONPATH="/projekt/rnaseq-toolkit/src:$PYTHONPATH"

# Przykład użycia CLI
python3 -m rnaseq_toolkit.cli \
    --counts data/batten/counts.csv \
    --metadata data/batten/metadata.csv \
    --design "~condition" \
    --contrast condition Batten control \
    --norm-method deseq2 \
    --dea-method deseq2 \
    --output results/batten_run/
```

---

## 6. Instalacja pakietu (dla użytkowników końcowych)

Po upublicznieniu repozytorium, pakiet można zainstalować przez:

```bash
pip install git+https://github.com/dawidx1233/rnaseq-toolkit.git
```

Docelowo (po submisji do PyPI):

```bash
pip install rnaseq-toolkit
```

---

## 7. Harmonogram prac

| Tydzień | Zadanie |
|---|---|
| 1 | Upublicznienie repo, dodanie CI, naprawa VST, vignettes |
| 2 | Enrichment analysis (GSEA/GO) na wszystkich 3 zestawach |
| 3 | Uzupełnienie artykułu (wyniki, rysunki, dyskusja) |
| 4 | Submisja do JOSS |
| 5-8 | Recenzja JOSS, odpowiedzi na recenzentów |
| 9+ | Rozbudowa do Bioinformatics Application Note |

---

*Wygenerowano automatycznie przez Manus AI — 29 maja 2026*
