# General
The following block has three `call{}` statements. In the first two we use the experimental data from GLDS to find upregulated genes and proteins. In the third we use these in conjucnction with the SPOKE reference network to map out the actual steps of mRNA transcription and protein translation and their connection to known diseases.

## Defining the database
In each call statement `USE` defines the database to search. `compositenasa.glds` holds the imported datasets from glds. `compositenasa.human` holds the a complete version of the SPOKE network.

## Collecting data  
Then we use `MATCH` to do several tasks at once, starting with identifying a study and any upregulated genes. The only information that is kept from this is the Gene node, by adding a variable, g, before the :Gene.  

Then the `WITH COLLECT` statement collects gene identifiers from the meta data for each gene node and stores it as glds_gids for use in the next call statement.  

In the second call statement we use the same methods but this time collect the upregulated protein identifiers.

In the third call statement we carry over the important gene information from before, but now point our search to the SPOKE reference database.

The `MATCH` function here describes a specific relatioship that connects our upregulated genes and proteins to diseases. 

This final matching expression requires that a gene encodes a protein. To satisfy this the gene and protein would've had to been upregulated at the transcrip and protein level. However the relationship connecting to the disease node is left blank `-[]->` which allows any of the connections between protein and disease in the SPOKE database to be returned. 


```Cypher
CALL {
USE compositenasa.glds
MATCH (n0:Study)-[]-(:Result)-[:UPREGULATION_RuMG]-(:MouseGene)-[]-(g:Gene)
WITH COLLECT(g.identifier) AS glds_gids
RETURN glds_gids
}

CALL{
USE compositenasa.glds 
MATCH (n0:Study {identifier:'GLDS-99'})-[]-(:Result)-[:UPREGULATION_RuMP]-(:MouseProtein)-[]-(p:Protein)
WITH COLLECT(p.identifier)AS glds_pids
RETURN glds_pids
}

UNWIND glds_pids AS glds_pid
UNWIND glds_gids AS glds_gid

CALL{
WITH glds_pid, glds_gid
USE compositenasa.human
MATCH path=(:Gene{identifier:glds_gid})-[ENCODES_GeP]->(:Protein{identifier:glds_pid})-[]->(:Disease)
RETURN path
}
RETURN path LIMIT 75
```