In this exmaple an in memory graph is projected from just the GLDS database searching for differentially regulated genes or proteins. Here we constrain the relationships but still allow flexibility. Additionally the log2fc property is collected from the relationships and stored as a variable (u) for later analysis.

```Cypher
USE compositenasa.glds
MATCH (r:Result)-[u:UPREGULATION_RuMP|DOWNREGULATION_RdMP|UPREGULATION_RuMG|DOWNREGULATION_RdMG]-(g:MouseProtein|MouseGene)
WITH gds.graph.project('gldsProteins8',r,g,
{
sourceNodeLabels:labels(r),
targetNodeLabels:labels(g),
relationshipProperties:u{.log2fc}
},
{undirectedRelationshipTypes: ['*']}
)AS a

RETURN a.graphName AS graph, a.nodeCount AS nodes, a.relationshipCount AS rels

```
image2

Above, we can see the relationship from study to dataset which in this case are different tissues. The different tissue then relate to their own set of differentially expressed transcripts, however some are shared and some are unique. Given the constraint of display all possible connections this image only represents a handful of connections. To understand which transcripts or proteins were shared among all datasets and studies we can utilze the graph functions built into Neo4j.

This simple example demonstrates the application of scoring nodes based on their calculated "betweeness centrality", a common knowledge graph statistic. Here we call an in memory graph, calculate the betweeness, and use `.mutate` to write the new data to nodes in the memory graph. Additionally the log2fc was used to weight in the interactions.
```Cypher
USE compositenasa.glds
CALL gds.betweenness.mutate('gldsProteins8',{ relationshipWeightProperty: 'log2fc',mutateProperty:'betweenness' })
YIELD centralityDistribution, nodePropertiesWritten
RETURN centralityDistribution.min AS minimumScore, centralityDistribution.mean AS meanScore, nodePropertiesWritten
```
image3

The results return transcripts ranked by betweeness score. In this case we can think of it as which transcripts are most connected and how important were they. Below is the actual representation of the top transcript and the number of tissues to which it was linked.

image4