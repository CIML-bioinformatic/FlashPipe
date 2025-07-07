---
output: html_document
---

# Documentation Développeur — FlashPipe

## Présentation générale

**FlashPipe** est un pipeline conçu pour le traitement et l’analyse de données générées par la technologie FlashFB5P-seq. Il est structuré autour de deux grandes phases :

* **Phase Python** : préparation de la structure de projet, organisation des fichiers d’entrée, vérification des paramètres fournis dans le fichier de configuration YAML.
* **Phase R** : traitement des résultats zUMIs et TRUST4, contrôle qualité (QC) des données, génération de visualisations et compilation d’un rapport HTML interactif (HTML).

Le tout est orchestré par **Snakemake**, assurant une exécution automatisée, traçable et reproductible.

Tous les scripts nécessaires au pipeline sont regroupés dans : `PROJECT_NAME/EXPERIENCE_NAME/03_Script/01_FlashPipe/`

---

## Architecture général du pipeline

| Étape               | Langage       | Rôle principal                                                        |
| ------------------- | ------------- | --------------------------------------------------------------------- |
| 1. Structure        | Python        | Création de l’arborescence, copie des fichiers, liens symboliques     |
| 2. zUMIs            | Snakemake     | Exécution du pipeline zUMIs (expression ARN)                          |
| 2. TRUST4           | Snakemake     | Exécution de l'outil pour l'analyse des récepteurs BCR/TCR            |
| 3. Contrôle Qualité | R             | Analyse des sorties zUMIs/TRUST4, nettoyage, tries, visualisations    |
| 4. Rapport final    | R             | Génération du rapport HTML avec des figures intégrées                 |

---

## Phase 1 — Préparation (Python)

### Script : `create_folder_structure.py`

**Emplacement** : `00_organizeStructure/`

Ce script crée et **met en place** automatiquement la **structure** complète du **projet** et **configure automatiquement** les **fichiers** nécessaires aux prochaines étapes du pipeline.

**Fonctionnalités clés** :

* Lecture et vérification des **paramètres** fournis dans le fichier **`config_FlashPipe.yml`** (type, format, logique).
* Génération d’un dictionnaire de variables pour copier.yml présent dans le répertoire du template.
* Personnalisation et **exécution du template** avec copier.yml présent dans le répertoire du template, pour attribuer les valeurs personnalisées au fichier.
* **Création** de **liens symboliques** vers les fichiers FASTQ.
* **Copie** des fichiers d’**IndexSort** et du **fichier GSF** dans les répertoires adéquats.

**Fichier associé** :

* `create_folder_structure_function.py` : contient l'ensemble des **fonctions** créées et utilisées dans le fichier create_folder_structure.py pour la **mise en place de la structure**, des **vérifications** de certains fichiers et des paramètres.

---

## Phase 2 — Contrôle Qualité (R)

Tous les scripts R employés pour le contrôle qualité sont situés dans `03_QC/`. Voici un aperçu synthétique de l'ensemble des fichiers et de leurs utilités :

| Fichier                      | Objectif                                                                                                            |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `00_generalDeps.R`           | Chargement des packages et définitions des fonctions communes utilisés (graphiques, traitement, etc.)               |
| `01_prepareData.R`           | Chargement de toutes les données d’entrée (zUMIs, TRUST4, IndexSort, GSF) pour vérifications (type, format...)      |
| `02_formatCountTable.R`      | Formatage des matrices de comptage, tri des données et génération de l'objet Seurat                                 |
| `03_computeERCCUMIPercent.R` | Calcul du pourcentage d’UMI assignés aux ERCC                                                                       |
| `04_computeERCCAccuracy.R`   | Évaluation de la corrélation entre les ERCC observés et attendus (vérifier l'étape d'amplification)                 |
| `05_computeRNAUMI.R`         | Analyse de la distribution des UMI ARN par puits                                                                    |
| `06_computeFeatureNumber.R`  | Analyse des nombres de gènes exprimés (features) pour toutes les plaques                                             |
| `07_computeIndexSort.R`      | Analyse des données de tri FACS (IndexSort), visualisation plaque par plaque                                        |
| `08_displayMetaData.R`       | Traitement et visualisation des métadonnées issues des fichiers GSF                                                 |
| `09_computeTrust4BCR.R`      | Analyse des isotypes BCR (classes et chaînes productives/non productives)                                           |
| `10_computeTrust4TCR.R`      | Analyse des isotypes TCR (alpha, beta, gamma, delta) et productivité des chaînes                                    |
| `analysisParams.R`           | Initialisation des constantes globales, chargement des options du YAML, palettes de couleur, noms de colonnes, etc. |

---

## Synthèse des Entrées et Sorties attendus

### Données d'entrées et leur positions

Répertoire : `EXPERIENCE_NAME/`

| Répertoire                    | Contenu attendu                   | Format      |
| ----------------------------- | --------------------------------- | ----------- |
| `00_RawData/00_RNA/`          | Fichiers FASTQ                    | `.fastq.gz` |
| `00_RawData/01_IndexSort/`    | Données FACS (IndexSort)          | `.csv`      |
| `01_Reference/00_Experiment/` | Plan de plaque (Metadata -- GSF)  | `.xlsx`     |
| `01_Reference/01_zUMIs/`      | Fichier de configuration zUMIs    | `.yml`      |
| `02_Container/`               | Conteneurs Singularity            | `.sif`      |
| `03_Script/01_FlashPipe/...`  | Scripts Python/R du pipeline      | `.py`, `.R` |
| `04_Workflow/01_snakemake/`   | Fichier Snakemake des deux étapes | `.yml`      |

---

### Fichiers produits

Répertoire : `EXPERIENCE_NAME/05_Output/01_FlashPipe/03_QC/`

| Étape/Script             | Type de données                                  | Format         |
| ------------------------ | ------------------------------------------------ | -------------- |
| `zUMIs`                  | Comptages UMI bruts                              | `.rds`         |
| `TRUST4`                 | Résultats BCR/TCR                                | `.tsv`         |
| `02_formatCountTable.R`  | Comptages RNA/ERCC, objets Seurat                | `.csv`, `.rds` |
| Scripts QC (`03` à `10`) | Figures descriptive (globaux et par plaque)      | `.png`, `.pdf` |
| Rapport HTML final       | Rapport interactif (figures, résumés)            | `.html`        |

---

### Figures générées

Graphique générés par les scripts présent dans le rapport HTML :

| Script                       | Figures construites                                      |
| ---------------------------- | -------------------------------------------------------- |
| `03_computeERCCUMIPercent.R` | Violin plot, plate plot                                  |
| `04_computeERCCAccuracy.R`   | Violin plot, plate plot, scatter plot                    |
| `05_computeRNAUMI.R`         | Violin plot, plate plot                                  |
| `06_computeFeatureNumber.R`  | Violin plot, UpSet plot, plate plot                      |
| `07_computeIndexSort.R`      | Violin plot ou Stacked Barchart (selon type de données)  |
| `08_displayMetaData.R`       | Violin plot ou Stacked Barchart (selon type de données)  |
| `09_computeTrust4BCR.R`      | Stacked Barchart, plate plot                             |
| `10_computeTrust4TCR.R`      | Stacked Barchart, plate plot                             |

---

#### Information complémentaire

* Les chemins sont construits dynamiquement depuis les fichiers `projectParams.R`, `experimentParams.R`, `analysisParams.R`, avec l'aide de Python, et du fichier de configuration (`config_FlashPipe.yml`).
* L’ensemble des détails des visualisations repose sur des constantes définis dans `00_generalDeps.R` (palettes de couleur, formats...).
* Les objets globaux (notamment les noms de colonnes, les dataframes triés, les noms des plaques...) sont stockés dans des listes et/ou des constantes durant l'analyse. Leur nommage et leur structure doivent rester cohérents pour être exploitables durant l'analyse.

---

