{{toc depth="1"/}}

La procédure d'analyse de données HiC (High Chromosome Contact) //in situ// présentée ici a pour objectif l'amélioration d'un assemblage génomique candidat à l'aide d'une carte HiC. Celle-ci représente les fréquences de contacts existants entre différentes régions de la chromatine, rendant compte de leur proximité physique. Les pipelines et outils utilisés dans cette procédure sont Juicer, 3D-DNA et Juicebox Assembly Tools développés par l' [[Aiden Lab>>https://github.com/aidenlab/]].

Une analyse de données HiC visant l'assemblage //de novo// d'un génome est également possible avec cette même suite d'outils.

Ces analyses nécessitent les ressources de calculs haute performance disponibles au Mésocentre Clermont Auvergne (HPC2).
Rappel pour demander un compte : [[Documentation Mésocentre Clermont Auvergne>>doc:Support informatique.Mésocentre Clermont-Auvergne||title="Documentation Mésocentre Clermont Auvergne"]]

= **Pipeline d'analyse de données HiC** =


(% style="color:blue" %)
== Prérequis ==

* **Juicer**
Copier le script permettant de générer //in silico//, pour l'assemblage candidat, les positions des sites de restriction correspondants aux enzymes utilisées pour produire les librairies HiC: [[generate_site_positions.py>>https://github.com/aidenlab/juicer/tree/master/misc]]
Copier le code du pipeline Juicer, dossier "SLURM -scripts for running pipeline and postprocessing on SLURM": [[https://github.com/aidenlab/juicer||title="https://github.com/aidenlab/juicer"]]
Paramétrer le script wrapper juicer.sh spécifiquement au serveur HPC2.

(((
{{code}}
## Juicer version 1.6
## Set the following variables to work with your system
load_bwa="module load bwa/0.7.17 java/oracle-1.8.0_45"
load_samtools="module load gcc/8.1.0 samtools/1.9"
load_gpu="module load cuda/10.2.89"
juiceDir="/storage/scratch/"$USER"/juicer-1.6"
## default queue, can also be set in options
queue="smp"
queue_time="3-00:00:00"
## default long queue, can also be set in options
long_queue="smp"
long_queue_time="3-00:00:00"
{{/code}}
)))

* **3D-DNA**
Copier le code du pipeline 3D-DNA permettant, dans cette procédure, de visualiser l'assemblage candidat en produisant la carte de contact HiC: [[https://github.com/aidenlab/3d-dna/||title="https://github.com/aidenlab/3d-dna/"]]

* **Juicebox**
Télécharger en local [[Juicebox Assembly Tools>>https://github.com/aidenlab/Juicebox/wiki/Download]]

* **Autres outils requis**
Java (version >= 1.7) sous Linux et Windows, Burrows-Wheeler Aligner (BWA)

(% style="color:blue" %)
== __Etape 1:__ Juicer ==

Juicer transforme les séquences HiC brutes en liste de contacts HiC (paires de positions génomiques qui sont adjacentes l'une de l'autre dans un espace en 3D au cours de l'expérimentation). Pour cela, les paires de séquences HiC sont alignées sur l'assemblage candidat, les lectures dupliquées sont supprimées, et les paires de lectures s'alignant plus de trois fois sont écartées. Les contacts Hi-C sont listés dans le fichier de sortie merged_nodups.txt.

NB: Pour permettre des modifications de l'assemblage à l'échelle des scaffolds lors de la 3ème étape de la procédure avec Juicebox Assembly Tools, et d'en assurer également la traçabilité, il convient de fournir à Juicer non pas le multifasta des pseudomolécules, mais un multifasta des scaffolds dans l'orientation et l'ordre suivi pour la construction des pseudomolécules.


==== **Index BWA** ====

Charger les modules requis et préparer le fasta index et l'index BWA.

(((
{{code}}
#!/bin/bash
module load gcc/4.8.4 bwa/0.7.17 samtools/1.3

samtools index draft.fa
cut -f1,2 /storage/scratch/$USER/juicer-1.6/references/draft.fa.fai > /storage/scratch/$USER/juicer-1.6/references/draft.fa.sizes

sbatch -p smp -c 1 --mem=32G --wrap="ml bwa; bwa index /storage/scratch/$USER/juicer-1.6/references/draft.fa"
{{/code}}
)))


==== **Générer les positions des sites de restriction** ====

Générer //in silico// les positions des sites de restriction, ici exemple pour le coktail d'enzymes Arima. Une ligne dans le fichier de sortie produit correspond à une séquence (fragment de restriction) dans l'assemblage candidat.

(((
{{code}}
/scripts/generate_site_positions.py Arima draft.fa /storage/scratch/$USER/juicer-1.6/references/draft.fa
{{/code}}
)))


==== **Repertoire d'analyse** ====

Mettre en place le répertoire de travail **juiceDir="/storage/scratch/"$USER"/juicer-1.6"** contenant:

* un dossier **/CANDIDATE_ASSEMBLY/fastq** contenant les séquences HiC brutes (convention de nommage _*R1*.fastq and _*R2*.fastq.gz).
* un dossier **/scripts** contenant tous les scripts du pipeline Juicer
* un dossier **/references** contenant le fasta de l'assemblage candidat, son index BWA, le fichier indiquant la taille des pseudomolécules (draft.fa.sizes), et le fichier draft_Arima.txt contenant les positions des sites de restriction.

==== **Lancer Juicer** ====

Il s'agit maintenant de lancer le script wrapper juicer.sh grâce à un script sbatch. Les jobs du pipeline Juicer seront ensuite soumis sur le serveur de calculs avec les ressources requises (dont, pour l'exemple ci-dessous, 16 CPU pour les jobs concernant l'alignement BWA).

Voici la description de l'utilisation du script juicer.sh (description complète: ./scripts/juicer.sh -h)

(((
{{code}}
Usage: juicer.sh [-g genomeID] [-d topDir] [-q queue] [-l long queue] [-s site] [-S stage] [-p chrom.sizes path] [-y restriction site file] [-z reference genome file] [-C chunk size] [-D Juicer scripts directory] [-Q queue time limit] [-L long queue time limit] [-b ligation] [-t BWA threads] [-A account name] [-e] [-h] [-f] [-j]
{{/code}}
)))

Editer un script sbatch **exemple_juicer_launch.sh** permettant de lancer le pipeline Juicer

(((
{{code}}
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --job-name=juicerCandidate
#SBATCH --partition=smp

#Set up the temporary directory
SCRATCHDIR='/storage/scratch/$USER/juicer-1.6'
mkdir -p -m 770 $SCRATCHDIR

cd $SCRATCHDIR

./scripts/juicer.sh \
-g "draft" \
-D $SCRATCHDIR \
-d $SCRATCHDIR/CANDIDATE_ASSEMBLY \
-s "Arima" \
-p $SCRATCHDIR/references/draft.fa.sizes \
-y $SCRATCHDIR/references/draft_Arima.txt \
-z $SCRATCHDIR/references/draft.fa \
-t 16 \
-e
{{/code}}
)))

Soumission du script sbatch **exemple_juicer_launch.sh**

(((
{{code}}
sbatch exemple_juicer_launch.sh
{{/code}}
)))


(% style="color:blue" %)
== __Etape 2:__ 3D-DNA ==

Lorsque l'on dispose d'un assemblage candidat à l'état de (quasi) pseudomolécules, un script du pipeline 3D-DNA permet de visualiser directement la carte des contacts HiC à partir de la liste des alignements dédupliqués et filtrés préalablement par Juicer, i.e. le fichier merged_nodups.txt.

Obtenir le fichier .assembly selon l'exemple:

(((
{{code}}
sbatch -p smp -c 1 --wrap="awk -f /storage/scratch/$USER/3d-dna/utils/generate-assembly-file-from-fasta.awk $SCRATCHDIR/references/draft.fa > $SCRATCHDIR/references/draft.assembly"
{{/code}}
)))

Obtenir la carte de contacts HiC (fichier .hic) selon l'exemple:

(((
{{code}}
sbatch -p smp -c 256 --mem=500G --wrap="ml java/oracle-1.8.0_45; export IBM_JAVA_OPTIONS='-Xmx500g -Xgcthreads256';
export _JAVA_OPTIONS='-Xmx500g -Xms500g';
/storage/scratch/$USER/3d-dna/visualize/run-assembly-visualizer.sh $SCRATCHDIR/references/draft.assembly $SCRATCHDIR/CANDIDATE_ASSEMBLY/aligned/merged_nodups.txt"
{{/code}}
)))


(% style="color:blue" %)
== __Etape 3:__ Juicebox Assembly Tools (JBAT) ==

Ouvrir l'application JBAT en local. Importer le fichier .hic (menu "Open", puis "File", "local"), puis le fichier .assembly (menu "Assembly", puis "Import Map Assembly").
Après avoir consulté ce tutoriel [[vidéo>>https://www.youtube.com/watch?v=Nj7RhQZHM18]], apporter les corrections nécessaires à la carte HiC directement dans l'application JBAT.

Le but est d'obtenir des fréquences de contacts homogènes formant une diagonale pour chaque chromosome.

Enfin, exporter le fichier .assembly intégrant ces modifications (menu "Assembly", puis "Export Assembly").


(% style="color:blue" %)
== __Etape 4:__ générer la nouvelle séquence ==

Il s'agit maintenant d'obtenir le fichier fasta du génome dont on vient de réviser l'assemblage, i.e. les séquences des pseudomolécules intégrant les modifications d'assemblage apportées sur la base de la carte HiC.

Le pipeline 3D-DNA propose la procédure suivante:

(((
{{code}}
./3d-dna/run-asm-pipeline-post-review.sh –r draft_reviewed.assembly draft.fa merged_nodups.txt
{{/code}}
)))

Une autre possibilité est de construire les nouvelles pseudomolécules en intégrant les informations du fichier draft_reviewed.assembly conjointement aux positions des streches de N présents dans la séquences draft.fa (script interne au GDEC nCount.pl). Cette alternative permet, lors de la fragmentation d'un scaffold avec l'outil JBAT, de faire coincider la coupure avec le début ou la fin d'un strech de N, de manière à assurer la conservation de l'intégrité des gènes. Suite au calcul des nouvelles coordonnées des fragments de scaffolds, les scripts internes au GDEC subseq.pl et buildPseudomol.pl permettent d'obtenir les séquences des pseudomolécules intégrant toutes les modifications d'assemblage.


(% style="color:blue" %)
== __Etape 5:__ validation ==

Afin de valider la nouvelle version d'assemblage ainsi obtenue, une carte HiC pourra être produite et visualisée avec JBAT en répétant les étapes 1 à 3.