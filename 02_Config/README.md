# README — Pipeline d’analyse FlashFB5p-seq

**Auteur :** Drystan THYARION 
**Contact :** drystan.thyarion@gmail.com 

---

## Introduction

Ce pipeline a été conçu pour automatiser l’analyse bioinformatique de données générées par la technique FlashFB5p-seq, une méthode de transcriptomique 5’ à l’échelle cellulaire unique ou en mini-bulk, développée au CIML. Cette technologie permet l’analyse conjointe du transcriptome et des répertoires BCR/TCR à partir de cellules triées par FACS.

Ce pipeline a pour objectif de répondre aux exigences de reproductibilité, modularité et traçabilité des analyses bioinformatiques, en s’appuyant sur des principes FAIR et une logique de science ouverte. Il génère à la fin un rapport HTML interactif de contrôle qualité intégrant données RNA-seq, FACS et répertoires immunitaires.

---

## Structure du projet

La structure est générée automatiquement à partir d’un template grâce à l’outil [`Copier`](https://copier.readthedocs.io/en/stable/). Cette section décrit la structure initiale (`01_Template`) et celle générée après exécution complète du pipeline.

### Template de base (`01_Template/`)

```
01_Template/
│
├── copier.yml
├── README.md
└── {{experience_name}}
    ├── 02_Container
    │   └── Contient les fichier .sif et nécéssaire pour les outils
    ├── 03_Script
    │	└── 01_FlashPipe
    │		└── Fichier de script employé pour la strucutre, et l'analyse QC
    ├── {{ output_name}}
    │	└── 01_FlashPipe
    │		├── 01_zUMIs
    │		│	└── {% yield plate_name from plate_names.split(',') %}{{ plate_name }}{% endyield %}
    │		├── 02_trust4
    │		├── 03_QC
    │		│	└── {% yield plate_name from plate_names.split(',') %}{{ plate_name }}{% endyield %}
    │		└── 04_Analysis
    ├── {{rawdata_name}}
    │	├── 00_RNA
    │	└── 01_IndexSort
    ├── {{workflow_name}}
    │	└── 01_snakemake
    │		├── config.yaml.jinja
    │		├── snakefile.yaml
    │		└── snakefile_copier.yaml
    └── {{reference_name}}
        ├── 00_Experiment
        │	├── 02_IMGT
        │	├── cell_barcode_well.csv
        │	└── ERCC_concentration.csv
        ├── 01_zUMIs
        │    └── {% yield plate_name from plate_names.split(',') %}{{ plate_name }}{% endyield %}
        │    	└── {{plate_name}}.yaml.jinja
        └── 02_trust4
```

**Ne pas modifier manuellement les fichiers générés dans le template.** Cela risquerait de casser l’exécution du pipeline.

---

### Dossier de configuration minimal (`02_Config/`)

```
02_Config/
│
└── PROJECT_NAME
    └── EXPERIMENT_NAME
        └── 01_Reference
            └── config_FlashPipe.yml
```

Ce dossier doit être copié pour chaque nouvelle analyse, avec adaptation du nom de projet, nom d’expérience, et remplissage du fichier `config_FlashPipe.yml`.

---

### Arborescence après exécution complète
(valable pour N plaques et toutes les options activées)

```
.../
│
└── PROJECT_NAME
    └── EXPERIMENT_NAME
    	├── 00_RawData
    	│	├── 00_RNA
    	│	│	└── Fichier symbolique des FastQ (Read 1 et 2) par plaque.
	│	└── 01_IndexSort
	│		└── 1 fichier d'index sort copier par plaque.
    	├── 01_Reference
    	│	├── 00_Experiment
    	│	├── 01_zUMIs
    	│	│	└── Plate_N
    	│	│		└── Plate_N.yaml
    	│	├── 02_trust4
    	│	└── config_FlashPipe.yml
    	├── 02_Container
    	├── 03_Script
	├── 04_Workflow
	│	├── config.yaml
    	│	├── snakefile.yaml
    	│	└── snakefile_copier.yaml
        └── 05_Output
        	├── 01_zUMIs
        	│	└── Plate_N
        	│		├── zUMIs_output
        	│		│	└── expression
        	│		│		├── Plate_N.dgecounts.rds
        	│		│		└── Autre fichier (.txt, .loom)
        	│		└── Fichier de sorties (.txt, .bam, .bai, .gtf, .tab ...)
        	├── 02_trust4
        	│	└── Plate_N
        	│		├── Plate_N_barcode_airr.tsv
        	│		└── Autre fichier (.tsv, .out, .fq, .fa)
        	├── 03_QC
        	│	├── PROJECT_NAME_EXPERIMENT_NAME_03_QC.html
        	│	├── Multiples graphiques sorties du QC (regroupant des informations généraux entre les plaques)
        	│	└── Plate_N
        	│		├── Multiples fichiers .csv provenant des données zUMIs retravailler (modifier)
        	│		├── Multiples graphique de la plaque N
        	│		└── Scatter_Plot_by_Well
        	│			└── Ensemble des graphique de corrélation pour tout les puits présent dans la plaque
        	└── 04_Analysis
```

- `Plate_N.yaml` un fichier de config pour l'outils zUMIs (il s'agit d'1 fichier par plaque)
- `Plate_N.dgecounts.rds` fichier rds de sortie de zUMIs qui sera utilisé dans l'analyse QC
- `Plate_N_barcode_airr.tsv` fichier de sortie de Trust4 qui seras utilisé pour la partie BCR et TCR de l'analyse de QC
- `PROJECT_NAME_EXPERIMENT_NAME_03_QC.html`  fichier/rapport de sortie de l'analyse de QC contenant l'ensemble des graphiques etc...

---

## Lancement du pipeline

### Étape 0 : Initialisation

1. Copier le dossier `02_Config` dans un nouveau projet.
2. Renommer `PROJECT_NAME` et `EXPERIMENT_NAME`.
3. Modifier et remplir `config_FlashPipe.yml` :
   - Noms des plaques
   - Mode (single-cell / mini-bulk)
   - Présence ou non de BCR/TCR
   - Chemins vers les fichiers FASTQ, FACS, GSF
   - Espèce (human / mouse)

### Étape 1 : Génération de la structure (Copier)

1. Se positionner dans EXPERIMENT_NAME.
2. Créer un environnement virtual :
   - virtualenv <NOM_ENVIRONNEMENT_VIRTUELLE>
   - source <NOM_ENVIRONNEMENT_VIRTUELLE>/bin/activate
   - pip install snakemake
   - pip install pulp==2.7.0
3. Lancer la commande : TEMPLATE_PATH=/mnt/DOSI/MASTER_TEMP/CB2M/Project/01_FlashPipe/01_Template
4. Lancer la commande : snakemake -j 1 --config template_path=${TEMPLATE_PATH} --snakefile "${TEMPLATE_PATH}/{{experience_name}}/{{workflow_name}}/01_snakemake/snakefile_copier.yaml" --use-singularity --singularity-args "-B /mnt:/mnt" --dryrun
5. Lancer la commande sans --dryrun, si aucune erreur n'est présente.

### Étape 2 : Utilisation des outils et générations de l'analyse QC
Lancer les outils et l'analyse QC :

1. Se connecter à un Gauss.
2. Se positionner dans EXPERIMENT_NAME.
3. Lancer la commande : snakemake -j 1 --snakefile 04_Workflow/01_snakemake/snakefile.yaml --use-singularity --singularity-args "-B /mnt:/mnt -B /tmp:/tmp" --dryrun
4. Lancer la commande sans --dryrun, si aucune erreur n'est présente.

# Puis sans --dryrun pour exécuter réellement



