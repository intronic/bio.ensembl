(defproject intronic/bio.ensembl "0.1.3-SNAPSHOT"
  :description "Clojure API for Ensembl data access"
  :url "intronic/bio.ensembl"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}

  :jvm-opts ^:replace ["-XX:MaxPermSize=256m" "-Xmx4g" "-server" "-d64" "-Djava.net.preferIPv4Stack=true"]

  :dependencies [[org.clojure/clojure "1.6.0"]
                 [org.clojure/algo.generic "0.1.2"]
                 [uk.ac.roslin/ensembl-data-access "1.17"]
                 [uk.ac.roslin/ensembl-config "1.71"]
                 [de.kotka/lazymap "3.1.1"]
                 [gavagai "0.3.1"]]

  :profiles {:dev {:source-paths ["dev"]
                   :dependencies [[org.clojure/tools.namespace "0.2.4"]]}}

  :source-paths ["src"]
  :test-paths ["test"]

  :lein-release {:deploy-via :clojars
                 :scm :git}

  :repositories { "jensembl" {:url "http://jensembl.sourceforge.net/m2-repo"
                              :checksum :ignore :snapshots false }
                  "biojava" {:url "http://www.biojava.org/download/maven/"
                             :checksum :ignore :snapshots false }
                  "snapshots" ~(str "file://" (System/getProperty "user.home") "/clj-repo/snapshots")
                  "releases" ~(str "file://" (System/getProperty "user.home") "/clj-repo/releases")}

  :warn-on-reflection false)
