(ns org.bioclojure.bio.ensembl.core-test
  (:require [clojure.test :as t :refer :all]
            [org.bioclojure.bio.ensembl.core :refer :all]))

(deftest test-coords
  (testing "round trip conversion"
    (is (= [10 20] (coord-vec (coord 10 20))))
    (is (= [10 20 1] (coord-vec (coord 10 20 1))))
    (is (= [10 20 -1] (coord-vec (coord 10 20 -1)))))
  (testing "strand"
    (is (thrown? clojure.lang.ExceptionInfo (strand+? nil)))
    (is (thrown? clojure.lang.ExceptionInfo (strand-? nil)))
    (is (= true (strand+? (nth (coord-vec (coord 10 20 1)) 2))))
    (is (= false (strand-? (nth (coord-vec (coord 10 20 1)) 2))))
    (is (= false (strand+? (nth (coord-vec (coord 10 20 -1)) 2))))
    (is (= true  (strand-? (nth (coord-vec (coord 10 20 -1)) 2)))))
  (testing "count"
    (is (= 1 (coord-count (coord 1 1))))
    (is (= 1 (coord-count (coord 1 1 1))))
    (is (= 1 (coord-count (coord 1 1 -1))))))
