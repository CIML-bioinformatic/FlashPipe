---
output: html_document
---

# Documentation pour la préparation des données expérimentales

Cette documentation détaille les étapes indispensables à suivre pour fournir les données brutes dans un format adéquat, en vue de l’analyse bioinformatique des données issues de **FlashFB5P-seq** par le logiciel **FlashPipe**.

---

## Étape 1 : Fichier GSF — Métadonnées

Les métadonnées doivent être fournies dans un fichier Excel au format **CSV**, appelé **GSF** (*General Submission Form*). Ce fichier doit être structuré selon les règles suivantes :

### Structure du fichier de la GSF

- Un onglet nommé `General Info` devant contenir, dans la section `Libraries information`, la liste des noms de plaques dans la colonne **Library ID**. Ces noms ne doivent contenir **ni caractères spéciaux**, **ni espaces** : utilisez uniquement des **lettres**, **chiffres**, le caractère **`_`** (underscore) et le caractère **`-`** (tiret).
- Un **onglet (feuille)** par **plaque**. Le **nom de chaque onglet** doit correspondre **exactement** au **nom de la plaque** qu’il représente comme indiqué dans l'onglet `General Info`.


### Structure d’un onglet (correspondant à une plaque)

Chaque onglet décrivant une plaque doit contenir un tableau représentant le **plan de la plaque**, avec les métadonnées associées à chaque puits.

#### Consignes :

- La **première colonne** doit obligatoirement s’appeler **`WellID`** (respecter les majuscules/minuscules). Dans cette colonne doit apparaitre l'ID du puits (A1, A2, A3...).
- Les **noms de colonnes** suivants ne doivent contenir **ni caractères spéciaux**, **ni espaces** : utilisez uniquement des **lettres**, **chiffres**, le caractère **`_`** (underscore).
- Le nom et l'ordre des colonnes est libre, mais **ces noms de colonnes doivent être les mêmes dans tous les onglets décrivant les plaques**.
- Les **cellules en dehors du tableau** (hors plan de plaque) doivent rester **vides** (il arrive parfois qu'en cas de copier/coller d'une GSF précédente, des cellules contiennent encore de l'information).
- Si **une donnée n’est disponible** pour un **puits**, laissez simplement la cellule **vide** (ne pas remplir avec des zéros ou des mentions du type `"NA"`). Si **un puits n'a aucune information associée**, vous pouvez soit créer une ligne vide en mettant l'ID du puits en premiere colonne et rien ensuite, soit ne pas faire apparaitre la ligne de ce puits.

#### Exemple du contenu d'un onglet :

```
| WellID | Sex | SortPheno | MouseID | Species | index | Age |
| ------ | --- | --------- | ------- | ------- | ----- | --- | 
| A1     | F   | indexsort | AER1658 | mouse   | ATTGC | 22  |
| A2     |     |           |         |         |       |     |
| A3     | F   | indexsort | AER1658 | mouse   | ATTGC | 22  |
| A4     | M   | indexsort | AER1769 | mouse   | AATCC | 23  | 
| A6     | M   | indexsort | AER1769 | mouse   | AATCG | 22  | 
| ...    | ... | ...       | ...     | ...     | ...   |...  |
```

Dans cet exemple :

- La **colonne WellID **est bien positionnée en **premier** dans le fichier, et **respecte la nomenclature**.
- Le **format est respecté** : noms de colonnes valides, structure homogène, aucune cellule superflue en dehors du tableau.
- les puits `A2` et `A5` n'ontpas de données. La ligne du puits `A2` est **laissé vide**, sans être supprimée du tableau alors que la ligne du puits `A5`n'apparait pas. Les deux choix sont possibles.


---

## Étape 2 : Fichier IndexSort — Données de FACS

Les données de FACS doivent être fournies dans de multiples fichiers Excel au format **CSV**. Ces fichiers doivent respecter les consignes suivante :

### Consignes :

- Les **données de FACS** doivent être **stockées** dans des **fichiers csv** pour **chaque plaque** (1 fichier csv par plaque).
- Les **fichiers** doivent porter la **nomenclature** suivante : `NOM_DE_PLAQUE`_indexsort.csv.

**Exemple** : une plaque portant le nom `P2_H9K0J2WY` aura une fichier nommé `P2_H9K0J2WY_indexsort.csv`

- La **première colonne** doit s’appeler **`WellID`** et contenir l'ID des puits  (A1, A2, A3...).
- Les nom des colonnes contenant des données **linéaires** (pas de transformation `asinh`) doivent commencer par **`Lin_`**, voir l'exemple ci-dessous
- Les **noms de colonnes** ne doivent contenir **ni caractères spéciaux**, **ni espaces** : utilisez uniquement des **lettres**, **chiffres**, le caractère **`_`** (underscore) et le caractère **`-`** (tiret).
- Le nom et l'ordre des colonnes est libre, **à condition que les noms soient strictement identiques entre les plaques** (donc entre les fichiers d'index sort).
- **Lister explicitement** les colonnes qui **ne sont pas** des données d’indexsort (par exemple : `Time`, `SortPheno`, etc.), pour les transmettre au bioinformaticien en charge de l'analyse.
- Les **cellules en dehors du tableau** (hors plan de plaque) doivent rester **vides**.
- Si **une donnée n’est disponible** pour un **puits**, laissez simplement la cellule **vide** (ne pas remplir avec des zéros ou des mentions du type `"NA"`). Si **un puits n'a aucune information associée**, vous pouvez soit créer une ligne vide en mettant l'ID du puits en premiere colonne et rien ensuite, soit ne pas faire apparaitre la ligne de ce puits.

#### Exemple :

Mon fichier contenant les données de FACS se nomme : `GEX_library_8_indexsort.csv`.

```
| WellID | Lin_FCS-A | CD19 | CD4     | Time | SortPheno |
|--------|-----------|------|---------|------|-----------|
| A1     | 526       | 120  | 142     | 48   | GC_LZ     |
| A2     | 872       | 89   | 245     | 55   | GC_DZ     |
| ...    | ...       | ...  | ...     | ...  | ...       |
```

Dans cet exemple :

- Le **nom du fichier** respecte la nomenclature imposée.
- La **colonne WellID **est positionnée en **premier** dans le fichier, et **respecte la nomenclature**.
- La **colonnes** qui contient des **données linéaires**, respecte la nomenclature (`Lin_` en préfixe).
- Le **format est respecté** : noms de colonnes valides, structure homogène, aucune cellule superflue en dehors du tableau.

---

## Étape 3 : Fichiers FastQ

Les fichiers FastQ issus du séquençage doivent respecter une nomenclature simple mais rigoureuse, à l’image des fichiers précédents. Normalement, ce nom est généré par le séquenceur à partir des noms de plaques fournis dans la GSF mais il est bon de vérifier avant analyse.

### Nomenclature :

Dans le cas du Read 1 :
`NOM_DE_PLAQUE`[...chaine libre...]`_`R1`_`[...chaine libre...].fastq.gz

Dans le cas du Read 2 :
`NOM_DE_PLAQUE`[...chaine libre...]`_`R2`_`[...chaine libre...].fastq.gz

### Règles :

- `R1` et `R2` sont **obligatoires** pour identifier les lectures forward et reverse.
- La partie variable (`[...]`) peut contenir des identifiants mais ne doivent contenir **ni caractères spéciaux**, **ni espaces** : utilisez uniquement des **lettres**, **chiffres**, le caractère **`_`** (underscore) et le caractère **`-`** (tiret).
- Les fichiers doivent être **compressés** au format `.gz`.

#### Exemples :

- `P2_H9K0J2WY_R1.fastq.gz`
- `P2_H9K0J2WY_R2.fastq.gz`

Ici, le nom de la plaque et les underscores (`_`) sont correctement positionnés.

- `P1_H9C0U2GX_S9_R1_L001.fastq.gz`
- `P1_H9C0U2GX_S9_R2_L001.fastq.gz`

Cet exemple est également valide : bien que des informations supplémentaires soient présentes, elles sont placées aux bons endroits.

---

