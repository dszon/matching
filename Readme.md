# matching
Compute a matching on a graph by an approximation scheme approach

# Compilation
Compile with
`make`

Run with
`matching --in filename1 --out filename2 [-q]`

filename1 is a comma separated csv containing an edge list of the graph. Option [-q] sqitches to silent output (only percentage of progress). Input format:

`nodeA,nodeB,weight`  
`nodeA,nodeC,weight`  
`...`  

where nodeA and nodeB are the node-IDs and weight is a float describing the edge weight.

filename2 is the list of the edges in a 2/3-optimal matching. The algorithm tries to maximize the sum of the weights.

Implemented is the approximation scheme of [Robert Preis](https://link.springer.com/chapter/10.1007%2F3-540-49116-3_24) and the augmentation algorithm due to [Doratha Drake and Stefan Hougardy](https://link.springer.com/chapter/10.1007/978-3-540-45198-3_2), which gives a linear time approximation scheme.
