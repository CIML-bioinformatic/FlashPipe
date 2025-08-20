# FlashFB5p-seq Pipeline – README

**Auteur** : Drystan THYARION
**Contact** : drystan.thyarion@gmail.com

---

## Objectif

Ce document décrit l’utilisation et la structure du pipeline d’analyse **FlashFB5p-seq**, basé sur **Snakemake** et utilisant des containers **Singularity**. Il est destiné à faciliter la mise en place, l'exécution et la compréhension du pipeline pour le traitement des données sortant du FlashFB5p-seq.

---

## Lancement du pipeline

### Étape 0 : Initialisation

1. Copier le __contenue__ du dossier `02_Config` dans un nouveau projet (et ne rien créer/ajouter de plus).
2. Renommer `EXPERIMENT_NAME` qui doit se trouver dans votre projet.
3. Modifier et remplir `config_FlashPipe.yml` :
   - Noms des plaques
   - Mode (single-cell / mini-bulk)
   - Présence ou non de BCR/TCR
   - Chemins vers les fichiers FASTQ, FACS, GSF
   - Espèce (human / mouse)

### Étape 1 : Génération de la structure (Copier)

1. Se connecter à un Gauss (serveur de calcul).
2. Se positionner dans EXPERIMENT_NAME.
3. Créer un environnement virtuel :
   - virtualenv <NOM_ENVIRONNEMENT_VIRTUELLE>
   - source <NOM_ENVIRONNEMENT_VIRTUELLE>/bin/activate
   - pip install snakemake
   - pip install pulp==2.7.0
4. Lancer la commande : TEMPLATE_PATH=/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/01_COMMON_DATA/04_TOOLS/01_FlashPipe/01_Template/
5. Lancer la commande : snakemake -j 1 --config template_path=${TEMPLATE_PATH} --snakefile "${TEMPLATE_PATH}/{{experience_name}}/04_Workflow/01_snakemake/snakefile_copier.yaml" --use-singularity --singularity-args "-B /mnt:/mnt" --dryrun
6. Lancer la commande sans --dryrun, si aucune erreur n'est présente.

### Étape 2 : Utilisation des outils et générations de l'analyse QC
Lancer les outils et l'analyse QC :

1. __Pré requis__ : Être connecté à Gauss (serveur de calcul), être positionné dans EXPERIMENT_NAME et avoir l'environnement virtuel activé.
2. Lancer la commande : snakemake -j 1 --snakefile 04_Workflow/01_snakemake/snakefile.yaml --use-singularity --singularity-args "-B /mnt:/mnt -B /tmp:/tmp" --dryrun
3. Lancer la commande sans --dryrun, si aucune erreur n'est présente.

---

## Fonctionnement du pipeline d'analyse de la technique de FlashFB5p-seq.

### Structure du Projet
Le répertoire Template sert de base pour générer une nouvelle expérience via l'outil **Copier**. 
Ne modifiez **aucun fichier manuellement**, sous peine de ne plus être en capacité de reproduire les analyses passées ou futures.

Le template sera copié, auquel celui sera appliqué des modifications mais sinon d'origine il est structuré de la manière suivante :

```
01_Template/
├── copier.yml                       # Configuration pour Copier
├── README.md
└── {{experience_name}}/
    ├── 02_Container/                # Fichiers Singularity (.sif) nécessaires
    ├── 03_Script/
    │	└── 01_FlashPipe/            # Fichier des scripts pour la structure du projet et de l'analyse du QC
    ├── {{output_name}}/
    │	└── 01_FlashPipe/
    │      ├── 01_zUMIs/
    		└── {% yield plate_name from plate_names.split(',') %}{{ plate_name }}{% endyield %}
    │      ├── 02_trust4/                        # Résultats Trust4
    │      ├── 03_QC/
    		└── {% yield plate_name from plate_names.split(',') %}{{ plate_name }}{% endyield %}
    │      └── 04_Analysis/                     # Analyse finale
    ├── {{rawdata_name}}/
    │   ├── 00_RNA/                  # FASTQ symboliques
    │   └── 01_IndexSort/            # Données de FACS
    ├── {{workflow_name}}/
    │	└── 01_snakemake/
    │      ├── config.yaml.jinja        # Fichier config pour le snakefile.yaml
    │      ├── snakefile.yaml           # Pipeline principal
    │      └── snakefile_copier.yaml    # Génération de la structure
    └── {{reference_name}}/
        ├── 00_Experiment/
        │	├── 02_IMGT/				# Références Trust4 (BCR/TCR)
        │	├── cell_barcode_well.csv
        │	└── ERCC_concentration.csv                
        ├── 01_zUMIs/
        │    └── {% yield plate_name from plate_names.split(',') %}{{ plate_name }}{% endyield %}
        │    	└── {{plate_name}}.yaml.jinja
        └── 02_trust4/
```

#### À noter : 
L'ensemble des fichiers présents dans le template ne doivent pas être supprimer ou modifier manuellement !!!!! (Si c'est le cas, alors, toutes les analyses précédentes ou futur ne pourront plus être faites, ou refaites).
- `copier.yml` est le fichier qui sera modifié par le pipeline lui même pour permettre de mettre en place la structure.
- `{{ }}` : éléments dynamiques substitués par Copier selon `copier.yml`.
- `snakefile_copier.yaml` est le fichier snakefile_copier, est le premier snakemake qui permet de lancer et de produire la strucutre de l'analyse.
- `config.yaml.jinja` fichier qui sera modifier par Copier, pour attribué les infos du fichier de config utilisateur, nécéssaire pour l'analyse QC et le fichier `snakefile.yaml`
- `snakefile.yaml` est le fichier qui permet de lancer l'analyse, et les outils.
- `experience_name` et `reference_name` etc... définissent le chemin de base pour les répertoires.
- `plate_names` est une chaîne de caractères contenant les noms des répertoires (en boucle) à créer, séparés par des virgules.
- `{% yield plate_name from plate_names.split(',') %}{{ plate_name }}{% endyield %}` permet de parcourir la liste `plate_names` pour créer les répertoires contenant le nom de chacune des plaques.
- `02_IMGT` contient les informations nécéssaire à l'outils Trust 4 (ce sont les genomes) pour permettre d'obtenir des infos sur les BCR et TCR.
- `{{plate_name}}.yaml.jinja` est le fichier qui sera modifier pour chaque plaque fournis, permettant à zUMIs de lire ces fichiers pour lancer les outils.

Le répertoires 02_Config, contiens les répertoires nécéssaire et initiaux pour lancer le pipeline (il suffit de copier le répertoires, en remplaçant les noms et le fichier config pour pouvoir lancer le code) :

---

```
02_Config/
└── PROJECT_NAME
    └── EXPERIMENT_NAME
        └── 01_Reference
            └── config_FlashPipe.yml
```

Pour la deuxième strucutre, je ne vais pas ré ecrire en détail toute les infos car certaines sont cités dans la précédentes strucutre. Cependant, voici les informations obtenus après que le pipeline fut lancer et terminé.
Nous prendrons comme exemple ici N plaques (sachant que N représente un ensemble de plaque possible) :

---

##### Remarque 
Il est important de prendre en compte, que l'arborescence présenter est celle qui possède l'ensemble des options activent. Si certaines options ne sont pas activé alors des résultats seront manquant. -
- **Mode Single-cell** : Toutes les analyses (zUMIs, FACS, Trust4) sont actives.
- **Mode Mini Bulk** : Analyse FACS désactivée.
- Si BCR/TCR sont désactivés, les dossiers Trust4 seront vides (sans erreur).

---

### Structure Générée par le Pipeline

```
.../
PROJECT_NAME/
    └── EXPERIMENT_NAME/
    	├── 00_RawData/
    	│   ├── 00_RNA/                 # Liens symboliques vers les FASTQ (Read 1 et 2) pour toutes les plaques renseignées.
    	│   └── 01_IndexSort/           # Fichiers FACS par plaque
    	├── 01_Reference/
    	│   ├── 00_Experiment/
    	│   ├── 01_zUMIs/
    	│   │	└── Plate_N
    	│   │		└── Plate_N.yaml
    	│   ├── 02_trust4/
    	│   └── config_FlashPipe.yml
    	├── 02_Container/
    	├── 03_Script/
    	├── 04_Workflow/
    	│   ├── config.yaml            # Fichier modifier par Copier avec les éléments rentrés dans le fichier de config (config_FlashPipe.yml) de l'utilisateur
    	│   ├── snakefile.yaml
    	│   └── snakefile_copier.yaml
    	└── 05_Output/
        	├── 01_zUMIs/
        	│   └── Plate_N/
        	│       ├── zUMIs_output/
        	│	│	└── expression
        	│	│		├── Plate_N.dgecounts.rds
        	│	│		└── Autre fichier (*.txt / *.loom)
        	│       └── Autre fichier de sortie (*.txt / *.bam / *.bai / *.gtf / *.tab / ...)
        	├── 02_trust4/
        	│   └── Plate_N/
        	│       ├── Plate_N_barcode_airr.tsv
        	│       └── Autre fichier (*.tsv / *.fq / *.fa)
        	├── 03_QC/
        	│   ├── <PROJECT_NAME>_<EXPERIMENT_NAME>_03_QC.html
        	│   ├── Graphiques globaux de l'analyse QC (multi-plaques)
        	│   └── Plate_N/
        	│       ├── Multiples fichiers *.csv provenant de zUMIs retravailler et modifier (mis au propre)
        	│       ├── Multiples graphiques pour chacune des plaques (dans notre cas N)
        	│       └── Scatter_Plot_by_Well/
        	│		└── Graphique de corrélation de tous les puits d'une plaque
        	└── 04_Analysis/


```

- `Plate_N.yaml` un fichier de config pour l'outils zUMIs (il s'agit d'1 fichier par plaque)
- `Plate_N.dgecounts.rds` fichier rds de sortie de zUMIs qui sera utilisé dans l'analyse QC
- `Plate_N_barcode_airr.tsv` fichier de sortie de Trust4 qui seras utilisé pour la partie BCR et TCR de l'analyse de QC
- `<PROJECT_NAME>_<EXPERIMENT_NAME>_03_QC.html`  fichier/rapport de sortie de l'analyse de QC contenant l'ensemble des graphiques etc...

---


