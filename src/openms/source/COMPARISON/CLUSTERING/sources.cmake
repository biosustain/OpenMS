### the directory name
set(directory source/COMPARISON/CLUSTERING)

### list all filenames of the directory here
set(sources_list
AverageLinkage.cpp
ClusterAnalyzer.cpp
ClusterFunctor.cpp
ClusterHierarchical.cpp
CompleteLinkage.cpp
EuclideanSimilarity.cpp
MultiplexCluster.cpp
MultiplexGrid.cpp
MultiplexLocalClustering.cpp
SingleLinkage.cpp
SILACClustering.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\COMPARISON\\CLUSTERING" FILES ${sources})

