In the Cypher query basics we searched two databases and asked a question of our reference database (SPOKE). A feature of working in CQL is to generate in memory graphs which hold a collection of defined variables and relationships found in experimental data. We can then use overlaps from experimental nodes and reference nodes as an entry point to the information in the reference graph.

image5

Graph projection relies on selecting data from our experimental data first to transfer to the referecne. The first portion is similar to the calls made in Cypher query basics.

```Cypher
CALL {
USE compositenasa.glds
MATCH (n0:Study)-[]-(:Result)-[:UPREGULATION_RuMG]-(:MouseGene)-[]-(g:Gene)
WITH COLLECT(g.identifier) AS glds_gids
RETURN glds_gids
}  

CALL{
USE compositenasa.glds 
MATCH (n0:Study)-[]-(:Result)-[:UPREGULATION_RuMP]-(:MouseProtein)-[]-(p:Protein)
WITH COLLECT(p.identifier)AS glds_pids
RETURN glds_pids
}  

UNWIND glds_pids AS glds_pid
UNWIND glds_gids AS glds_gid
```
The third call block again uses the experimental gene and protein ID's but also creates a new variables to store the node information from the SPOKE network: g0, p, and d. Notably, this query does not specify a direction or type of relationship, which will return all possible connections.

Node information is then saved as a new in memory graph database "geneProtDisease". Here it is converted to a variable "a" and is `RETURN` inside this `CALL` function delivers the graph.
```Cypher
CALL{
WITH glds_pid, glds_gid
USE compositenasa.human

MATCH (g0:Gene{identifier:glds_gid})--(p:Protein{identifier:glds_pid})--(d:Disease)--(:Gene)

WITH gds.graph.project('geneProtDisease',g0,d,{sourceNodeLabels:labels(g0),targetNodeLabels:labels(d)})AS a
RETURN a
}
```
Essentially, this uses data from our experimental database and links to our reference database. By creating an in memory graph of the experimental database we can prune some of the more distant relationships that we are not interested in.

image6

Finally we can call the new graph and display the experimental data using the connections in the SPOKE database.
```Cypher
RETURN a.graphName AS graph, a.nodeCount AS nodes, a.relationshipCount AS rels
```

![graphProjection_example1](https://github.com/user-attachments/assets/243bd273-6794-4c5b-8eea-0247268f8b8a)
