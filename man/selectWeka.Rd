\name{selectWeka}
\alias{selectWeka}
\title{Feature Selection by Weka}
\description{
  feature selection by Weka.
}
\usage{
  selectWeka(train, evaluator="CfsSubsetEval", search="BestFirst", n)  
}
\arguments{
  \item{train}{a data frame including the feature matrix and class label of training set.}
  \item{evaluator}{a string for the feature selection method used by WEKA. This 
    must be one of the strings "CfsSubsetEval", "ChiSquaredAttributeEval", 
    "InfoGainAttributeEval", or "SVMAttributeEval".} 
  \item{search}{a string for the search method used by WEKA. This must be one 
    of the strings "BestFirst" or "Ranker".} 
  \item{n}{an integer for the number of selected features.}   
}
\details{  
  Parameter "evaluator" supportes three feature selection methods provided by WEKA:
  "CfsSubsetEval": Evaluate the worth of a subset of attributes by considering 
                   the individual predictive ability of each feature along with 
                   the degree of redundancy between them.
  "ChiSquaredAttributeEval": Evaluate the worth of an attribute by computing the 
                             value of the chi-squared statistic with respect to 
                             the class.
  "InfoGainAttributeEval": Evaluate attributes individually by measuring 
                           information gain with respect to the class.
  "SVMAttributeEval": Evaluate the worth of an attribute by using an SVM classifier. 
                      Attributes are ranked by the square of the weight assigned 
                      by the SVM. Attribute selection for multiclass problems is 
                      handled by ranking attributes for each class seperately 
                      using a one-vs-all method and then "dealing" from the top 
                      of each pile to give a final ranking.
                      
  Parameter "search" supportes three feature subset search methods provided by WEKA:
  "BestFirst":  Searches the space of attribute subsets by greedy hillclimbing 
                augmented with a backtracking facility. Setting the number of 
                consecutive non-improving nodes allowed controls the level of 
                backtracking done. Best first may start with the empty set of 
                attributes and search forward, or start with the full set of 
                attributes and search backward, or start at any point and search 
                in both directions (by considering all possible single attribute 
                additions and deletions at a given point).
  "Ranker": Ranks attributes by their individual evaluations.
 
}
\author{Hong Li}

  
