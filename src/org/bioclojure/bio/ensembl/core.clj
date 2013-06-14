(ns org.bioclojure.bio.ensembl.core
  (:require [clojure.java.io :as io]
            [clojure.string :refer (upper-case)]
            [org.bioclojure.bio.ensembl.config :only (data-source)])
  (:import [uk.ac.roslin.ensembl.config RegistryConfiguration DBConnection$DataSource]
           [uk.ac.roslin.ensembl.dao.database DBRegistry DBSpecies]
           [uk.ac.roslin.ensembl.model Coordinate StableID]
           [uk.ac.roslin.ensembl.model.database Registry]
           [uk.ac.roslin.ensembl.model.core
            DNASequence Feature Gene Species Transcript Translation Chromosome]
           [uk.ac.roslin.ensembl.model.variation Variation]
           [uk.ac.roslin.ensembl.datasourceaware.core DATranslation DAGene DATranscript]))

(defonce ^:dynamic ^Registry *registry* nil)

(defmacro with-registry
  [registry & body]
  `(binding [*registry* ~registry] ~@body))

(defn set-registry!
  [registry]
  (alter-var-root #'*registry*
                  (constantly registry)
                  (when (thread-bound? #'*registry*)
                    (set! *registry* registry))))

(defn local-config
  [conf-file]
  (DBRegistry. (doto (RegistryConfiguration.)
                 (.setDBByFile (io/file conf-file)))
               true))

(defn data-source
  [ds]
  (DBConnection$DataSource/valueOf
   (.toUpperCase (name ds))))

(defn registry
  [ds]
  (DBRegistry. (data-source ds)))

(defn ensembl-version
  "Return current Ensembl version of connected database registry for species-name."
  [species-name]
  (-> *registry* (.getDatabase (name species-name)) .getDBVersion ))


(defn- list-species-transform
  [style]
  (get {:binomial #(.getSpeciesBinomial ^Species %)
        :common   #(.getCommonName ^Species %)
        :compara  #(.getComparaName ^Species %)
        :database #(.getDatabaseStyleName ^Species %)
        :display  #(.getDisplayName ^Species %)
        :short    #(.getShortName ^Species %)}
       style
       identity))

(defn list-species
  [& [style]]
  (map (list-species-transform style) (.getSpecies *registry*)))

(defn species
  "Ensembl species from name or keyword"
  ^Species
  [species-name]
  (or (.getSpeciesByEnsemblName *registry* (name species-name))
      (.getSpeciesByAlias *registry* (name species-name))))

(defn species-version
  "Return genome assembly version of species."
  [species-name]
  (.getAssembly (species species-name)))

(defn list-chromosomes
  [species-name]
  (map #(.getChromosomeName ^Chromosome %) (vals (.getChromosomes (species species-name)))))

(defn chromosome
  ^Chromosome
  ([species-name chromosome-name]
     (.getChromosomeByName (species species-name) chromosome-name))
  ([species-name chromosome-name ens-version]
     (.getChromosomeByName (species species-name) chromosome-name (str ens-version))))

(defn chromosome-dna-str
  "Return DNA for chromosome location."
  ([^Chromosome chromosome start end]
     (-> chromosome (.getSequenceAsString (int start) (int end))))
  ([^Chromosome chromosome position upstream downstream]
     (-> chromosome (.getSequenceAsString (int (- position upstream))
                                          (int (+ position downstream))))))

(defn genes-on-region
  ([species-name chromosome-name begin end]
     (genes-on-region (chromosome species-name chromosome-name) begin end))
  ([^Chromosome chromosome begin end]
     (.getGenesOnRegion chromosome (int begin) (int end)))
  ([chromosome pos]
     (genes-on-region chromosome pos pos)))

(defn variations-on-region
  ([species-name chromosome-name begin end]
     (variations-on-region (chromosome species-name chromosome-name) begin end))
  ([^Chromosome chromosome begin end]
     (.getVariationsOnRegion chromosome (int begin) (int end)))
  ([chromosome pos]
     (variations-on-region chromosome pos pos)))

(defn gene
  ([species-name gene-stable-id]
     (.getGeneByStableID (species species-name) gene-stable-id))
  ([species-name gene-stable-id ens-version]
     (.getGeneByStableID (species species-name) gene-stable-id (str ens-version))))

(defn gene-transcripts
  ([^Gene gene]
     (.getTranscripts gene)))

(defn gene-canonical-transcript
  ([^Gene gene]
     (.getCanonicalTranscript gene)))

(defn transcript
  "Get transcript by stable ID"
  ([species-name transcript-stable-id]
     (.getTranscriptByStableID (species species-name) transcript-stable-id))
  ([species-name transcript-stable-id ens-version]
     (.getTranscriptByStableID (species species-name) transcript-stable-id (str ens-version))))

(defn transcript-exons
  "Get exons from transcript"
  ([species-name transcript-stable-id]
     (transcript-exons (transcript species-name transcript-stable-id)))
  ([transcript]
     (.getExons ^Transcript transcript)))

(defn gene-stable-id
  [^DAGene gene]
  (.getStableID gene))

(defn gene-name
  [^DAGene gene]
  (.getDisplayName gene))

(defn gene-description
  [^Gene gene]
  (.getDescription gene))

(defn transcript-stable-id
  [^DATranscript transcript]
  (.getStableID transcript))

(defn translation-transcript-stable-id
  [^DATranslation translation]
  (-> translation (.getTranscript) (.getStableID)))

(defn translation-stable-id
  [^DATranslation translation]
  (.getStableID translation))

(defn transcript-strand
  "Strand of transcript."
  [^Transcript transcript]
  (-> transcript (.getChromosomeMapping) (.getTargetCoordinates) (.getStrand)))

(defn coord-contains-point?
  "True if coord contains point"
  [^Coordinate coord pos]
  (.containsPoint coord (int pos)))

(defn coord-vec
  "Return a vector of coordinates [chr start end strand]."
  [^Coordinate coord]
  ((juxt #(.getStart ^Coordinate %) #(.getEnd ^Coordinate %) #(.getStrandInt ^Coordinate %)) coord))


(defn transcript-coord
  "Genomic coordinates of transcript start/stop/strand."
  [^Transcript transcript]
  (-> transcript (.getChromosomeMapping) (.getTargetCoordinates)))

(defn transcript-chrom
  "Genomic coordinates of transcript start/stop/strand."
  [^Transcript transcript]
  (-> transcript (.getChromosomeMapping) (.getTarget)))

(defn transcript-biotype
  [^DATranscript transcript]
  (.getBiotype transcript))

(defn transcript-gene
  [^Transcript transcript]
  (.getGene transcript))

(defn exon-rank
  "Rank for an exon"
  [exon]
  (and exon (-> exon .getRank)))

(defn exon-coord
  "Chromosome coordinates for an exon"
  ^Coordinate
  [exon]
  (-> exon .getChromosomeMapping .getTargetCoordinates))

(defn exon-coord-rank-vec
  "Chromosome coordinates and rank for an exon"
  [exon]
  (conj (coord-vec (exon-coord exon)) (exon-rank exon)))

(defn translation-transcript
  [^Translation translation]
  (.getTranscript translation))

(defn translation-gene
  [^Translation translation]
  (-> translation (translation-transcript) (transcript-gene)))

;;; Strand predicates
(defn strand-int+?
  [strand]
  (= 1 strand))

(defn strand-int-?
  [strand]
  (= -1 strand))

(defn strand+?
  [strand]
  (= uk.ac.roslin.ensembl.model.Coordinate$Strand/FORWARD_STRAND strand))

(defn strand-?
  [strand]
  (= uk.ac.roslin.ensembl.model.Coordinate$Strand/REVERSE_STRAND strand))

;;; Translation 
(defn transcript-canonical-translation
  "Returns the canonical translation for this transcript (if any)."
  [^Transcript transcript]
  (.getCanonicalTranslation transcript))

(defn transcript-translations
  "Returns the translations for this transcript (if any)."
  [^Transcript transcript]
  (.getTranslations transcript))

(defn protein-sequence
  "Returns the protein sequence for this translation."
  ^org.biojava3.core.sequence.ProteinSequence
  [^Translation translation]
  (.getProteinSequence translation))

(defn protein-aa
  "Returns the amino acid for this protein"
  [^org.biojava3.core.sequence.ProteinSequence protein aa-pos]
  (.getCompoundAt protein  aa-pos))

(defn protein-aa-str
  "Returns the amino acid for this protein"
  [^org.biojava3.core.sequence.ProteinSequence protein aa-pos]
  (str (.getCompoundAt protein  aa-pos)))

(defn str-protein-sequence
  "Returns a string of the protein sequence for this translation."
  ^String
  [^DATranslation translation]
  (.getProteinSequenceAsString translation))

(defn aa->chromosome
  "Convert AA position to chromosome location"
  [^DATranslation translation pos]
  (.getChromosomePositionFromAA translation (int pos)))

(defn aa<-chromosome
  "Convert AA position from chromosome location"
  [^DATranslation translation pos]
  (.getAAPositionFromChromosome translation (int pos)))

(defn cds<-chromosome
  "Get position relative to CDS from chromosome location."
  [^DATranslation translation pos]
  (.getBasePositionFromChromosome translation (int pos)))

(defn cds->chromosome
  "Get chromosome location from position relative to CDS."
  [^DATranslation translation pos]
  (.getChromosomePositionFromBASE translation (int pos)))

(defn transcript<-aa
  "Codon coordinates of amino acid relative to processed transcript."
  ^Coordinate
  [^DATranslation translation aa-pos]
  (.getProcessedTranscriptPositionFromAA translation (int aa-pos)))

(defn transcript-cds-start-coord
  "Start codon coordinate relative to processed transcript."
  ^Coordinate
  [translation]
  (transcript<-aa translation  1))

(defn transcript-cds-start
  "Start codon position relative to processed transcript."
  [^DATranslation translation]
  (.getStart (transcript-cds-start-coord translation)))

(defn cds-dna
  "Return BioJava DNASequence object for translated (CDS) region."
  ^DNASequence
  [^DATranslation translation]
  (.getTranslatedSequence translation))

(defn cds<-aa
  "Codon start position of amino acid relative to CDS start (1-based)."
  [aa-pos]
  (inc (* 3 (dec aa-pos))))

(defn aa-dna
  "Return DNA from CDS for amino acid position (1-based aa-pos)."
  [translation aa-pos]
  (let [pos (cds<-aa aa-pos)]
    (when-let [dna (cds-dna translation)]
      #_(assert (and (<= 1 aa-pos (.length (str translation)))
                   (pos? pos) (<= (+ pos 2) (.length (str dna))))
              (str "aa-pos=" aa-pos " length(aa)=" (.length (str translation))
                   " pos=" pos " length(dna)=" (.length (str dna))
                   (if (>= (.length (str translation)) 1)
                     (str " AA[" (.length (str translation))
                          "]=" (subs (str translation)
                                     (dec (.length (str translation))))))
                   (if (>= (.length (str dna)) 3)
                     (str " codon[" (/ (.length (str dna)) 3) "]=" (subs (str dna)
                                                                         (- (.length (str dna)) 3))))))
      (if (<= (+ pos 2) (.getLength dna))
        (-> dna (.getSubSequence (int pos) (int (+ pos 2))) (.getSequenceAsString))))))

(defn codon-dna
  "DNA of codon at chromosome position chrom-pos in translation."
  [translation chrom-pos]
  (aa-dna translation (aa<-chromosome translation chrom-pos)))

(defn cds-length
  "CDS length in processed transcript (including start/stop)."
  [^DATranslation translation]
  (.getLength translation))

(defn str-cds-dna
  "Return string of translated (CDS) region."
  ^String
  [translation]
  (.getSequenceAsString (cds-dna translation)))

(defn seq-cds-dna
  "Return seq of BioJava NucleotideCompound for translated (CDS) region."
  [translation]
  (seq (cds-dna translation)))
  
(defn aa-length
  "Amino acid length (including start/stop)."
  [translation]
  (/ (cds-length translation) 3))

(defn ccds-id
  [^Transcript transcript]
  (.getCcdsID transcript))

(defn highest-ensembl-version
  "Return the highest ensembl schema version"
  []
  (.getHighestEnsemblSchemaVersion *registry*))

(defn ensembl-versions
  "Return a seq of ensembl schema versions in the registry."
  []
  (seq (.getKnownSchemaVersions ^DBRegistry *registry*)))

(defn database
  "Get a single species core database version (or the most recent)
   Usage:
    (get-database \"human\")
    (get-database \"human\" 70)
    (get-database \"human\" \"70\")"
  ([species-name]
     (database species-name (highest-ensembl-version)))
  ([species-name ens-version]
     (.getDatabase *registry* species-name (str ens-version))))

(defn dna-complement
  "Complement of a org.biojava3.core.sequence.DNASequence"
  [^org.biojava3.core.sequence.DNASequence dna-seq]
  (-> dna-seq (.getComplement) (.getViewedSequence)))

(defn exon<-chromosome
  "Get exon from transcript and chromosome location"
  [transcript pos]
  (let [e (->> transcript
               transcript-exons
               (filter #(-> % exon-coord (coord-contains-point? pos))))]
    (assert (<= (count e) 1) (str "count=" (count e)))
    (first e)))

(comment

  (def ensreg (registry :ensembldb))

  (with-registry ensreg
    (list-species)
    (chromosome "human" "20")
    (gene "human" "ENSG00000153551"))

  (set-registry! ensreg)

  (list-species)

  (chromosome "human" "20")

  (list-chromosomes "human")

  (genes-on-region "human" "20" 1 100000)

  (genes-on-region (chromosome "human" "20") 1 100000)

  (ensembl-versions)

  
  (. (chromosome "human" "14") getSequenceAsString (int 81610500) (int 81610540))
  ;; "AAACGCCAGGCTCAGGCATACCGGGGGCAGAGGGTTCCTCC"
  
  ;; local Ensembl connections - for file format see:
  ;;   http://jensembl.svn.sourceforge.net/viewvc/jensembl/trunk/EnsemblTest/src/main/resources/
  (set-registry! (local-config "example_local_configuration.properties")))
