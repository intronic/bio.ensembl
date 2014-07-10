(defproject org.clojars.intronic/bio.ensembl "0.1.2-SNAPSHOT"
  :description "Clojure API for Ensembl data access"
  :url "org.clojars.intronic/bio.ensembl"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :repositories [["jensembl" {:url "http://jensembl.sourceforge.net/m2-repo" :checksum :ignore :snapshots false }]
                 ["biojava" {:url "http://www.biojava.org/download/maven/" :checksum :ignore :snapshots false }]]
  :dependencies [[org.clojure/clojure "1.4.0"]
                 [org.clojure/algo.generic "0.1.1"]
                 [uk.ac.roslin/ensembl-data-access "1.15"]
                 [uk.ac.roslin/ensembl-config "1.70"]
                 [de.kotka/lazymap "3.1.1"]
                 [gavagai "0.3.1"]
                 ]
  :warn-on-reflection true)
