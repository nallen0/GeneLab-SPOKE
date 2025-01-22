# SPOKE and NASA GeneLab Integration – Notes and Guide  

The SPOKE database, hosted by Neo4j, can be used as an analytical reference for operational datasets containing experimental data. The Neo4j platform is used to interface with the vast collection of connected SPOKE reference databases with the Cypher Query Language, a language based on SQL but specific for graph databases. Cypher is designed to use descriptive text to represent between nodes.
```
(node1) –[:relationships]-> (node2)
```
  Node types can be assigned containing different classes of data (e.g., Mission, Study, Results, Genes, Proteins) and each node defined by a label contains {properties}. Similarly, relationships (e.g., Up Regulation, Down Regulation) contain properties such as log2FC and Adjusted p value. 

  ### Examples
  1_Cypher_query_basics.md
  2_Graph_projection.md

  ### **Resources**
  See [Cypher Getting Started](https://neo4j.com/docs/getting-started/cypher-intro/) for an introduction to Cypher language or the [Cypher Manual](https://neo4j.com/docs/cypher-manual/current/introduction/) for complete documentation. 
  
  Visit the [Neo4j Graph Academy](https://graphacademy.neo4j.com/?_gl=1*lgv80f*_ga*MTQ4MjQ0Njg1NS4xNjkwMTUzNTQ3*_ga_DL38Q8KGQC*MTY5MTE3MDEzMy4yMi4xLjE2OTExNzIwMzUuNTguMC4w&_ga=2.106043587.1557279067.1691170135-1482446855.1690153547) for foundational courses.

  Additionally, the [Neo4j Sandbox](https://neo4j.com/sandbox/) provides guided tutorials working in projects querying with both in Cypher and a graphical user interface. This approach was the fastest and most intuitive way to learn Cypher and the capabilities of Neo4j, without reading the full manual. 

  ### Cypher and Graph Database Structures  

  The Cypher Query Language was designed for graph databases and to be easily adapted to for those with experience in SQL database queries. The graphical nature of the language and relationships makes it simple to learn even without database management experience. From a technical standpoint, Cypher is pitched as being a much simpler and faster form of exploring relationships in databases. The primary example from Neo4Jj is Cypher requiring significantly less lines of code to achieve the same relationship expression as an SQL query. This is achieved with nodes connected by edges, which themselves contain relational data and properties. Edges are the `[:relationships]` used to define queries while nodes contain specific data or results. 

  Highly connected data that is spread across many tables, such as in a multi-omics investigation, is well suited to graphical database representation. In such an investigation we are interested in end to end relationships of multiple datasets and need to quickly search for connections across tables. Graphical databases also enable analysis of such connections and even prediction of new ones based on existing relationships. Finally graph databases and concurrent analyses can be rapidly updated with additional data and relationships. 
  