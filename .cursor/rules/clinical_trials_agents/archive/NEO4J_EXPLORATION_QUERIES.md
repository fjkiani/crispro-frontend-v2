# 🔍 NEO4J EXPLORATION QUERIES

**Purpose:** Queries to explore the clinical trials graph database  
**Database:** `neo4j` (default)  
**Connection:** `neo4j+s://9669e5f3.databases.neo4j.io`

---

## **📊 OVERVIEW QUERIES**

### **1. Node Counts (Total Statistics)**
```cypher
MATCH (t:Trial) 
WITH count(t) AS trials
MATCH (o:Organization) 
WITH trials, count(o) AS orgs
MATCH (s:Site) 
WITH trials, orgs, count(s) AS sites
MATCH (p:PI) 
WITH trials, orgs, sites, count(p) AS pis
MATCH (c:Condition)
RETURN 
  trials AS `Total Trials`,
  orgs AS `Total Organizations`,
  sites AS `Total Sites`,
  pis AS `Total PIs`,
  count(c) AS `Total Conditions`
```

### **2. Relationship Counts**
```cypher
MATCH ()-[r:SPONSORS]->()
WITH count(r) AS sponsors
MATCH ()-[r2:CONDUCTED_AT]->()
WITH sponsors, count(r2) AS sites_rel
MATCH ()-[r3:LEADS]->()
WITH sponsors, sites_rel, count(r3) AS leads
MATCH ()-[r4:COLLABORATES_ON]->()
RETURN 
  sponsors AS `Sponsor Relationships`,
  sites_rel AS `Site Relationships`,
  leads AS `PI Lead Relationships`,
  count(r4) AS `Collaborator Relationships`
```

---

## **🔬 DETAILED EXPLORATION**

### **3. Sample Trials with Full Details**
```cypher
MATCH (t:Trial)
RETURN 
  t.nct_id AS NCT_ID,
  t.title AS Title,
  t.status AS Status,
  t.phase AS Phase
LIMIT 10
```

### **4. Trial with All Relationships (Single Trial Deep Dive)**
```cypher
MATCH (t:Trial {nct_id: 'NCT00000001'})
OPTIONAL MATCH (t)<-[:SPONSORS]-(sponsor:Organization)
OPTIONAL MATCH (t)<-[:COLLABORATES_ON]-(collab:Organization)
OPTIONAL MATCH (t)-[:CONDUCTED_AT]->(site:Site)
OPTIONAL MATCH (t)-[:LEADS]-(pi:PI)
OPTIONAL MATCH (t)-[:TREATS]->(condition:Condition)
RETURN 
  t.nct_id AS Trial_ID,
  t.title AS Title,
  collect(DISTINCT sponsor.name) AS Sponsors,
  collect(DISTINCT collab.name) AS Collaborators,
  collect(DISTINCT site.facility + ', ' + site.city + ', ' + site.state) AS Sites,
  collect(DISTINCT pi.name) AS Principal_Investigators,
  collect(DISTINCT condition.name) AS Conditions
LIMIT 1
```

### **5. All Organizations and Their Trial Counts**
```cypher
MATCH (o:Organization)-[:SPONSORS]->(t:Trial)
RETURN 
  o.name AS Organization,
  o.type AS Type,
  count(t) AS Trial_Count
ORDER BY Trial_Count DESC
LIMIT 20
```

### **6. Sites by State (Location Analysis)**
```cypher
MATCH (s:Site)
RETURN 
  s.state AS State,
  count(s) AS Site_Count,
  collect(DISTINCT s.city)[0..5] AS Sample_Cities
ORDER BY Site_Count DESC
LIMIT 20
```

### **7. Trials by Status**
```cypher
MATCH (t:Trial)
RETURN 
  t.status AS Status,
  count(t) AS Count
ORDER BY Count DESC
```

### **8. Trials by Phase**
```cypher
MATCH (t:Trial)
RETURN 
  t.phase AS Phase,
  count(t) AS Count
ORDER BY Count DESC
```

---

## **🔗 RELATIONSHIP EXPLORATION**

### **9. Organization → Trial → Site Chain**
```cypher
MATCH (org:Organization)-[:SPONSORS]->(t:Trial)-[:CONDUCTED_AT]->(site:Site)
RETURN 
  org.name AS Sponsor,
  t.nct_id AS Trial_ID,
  t.title AS Trial_Title,
  site.facility AS Site_Facility,
  site.city AS City,
  site.state AS State
LIMIT 20
```

### **10. Find Trials by Location (State-based)**
```cypher
MATCH (t:Trial)-[:CONDUCTED_AT]->(s:Site {state: 'CA'})
RETURN 
  t.nct_id AS NCT_ID,
  t.title AS Title,
  s.facility AS Facility,
  s.city AS City
LIMIT 20
```

### **11. Organizations Sponsoring Multiple Trials**
```cypher
MATCH (o:Organization)-[:SPONSORS]->(t:Trial)
WITH o, count(t) AS trial_count
WHERE trial_count > 1
RETURN 
  o.name AS Organization,
  o.type AS Type,
  trial_count AS Trials_Sponsored
ORDER BY trial_count DESC
```

### **12. Sites Hosting Multiple Trials**
```cypher
MATCH (s:Site)<-[:CONDUCTED_AT]-(t:Trial)
WITH s, count(t) AS trial_count
WHERE trial_count > 1
RETURN 
  s.facility AS Facility,
  s.city AS City,
  s.state AS State,
  trial_count AS Trials_Hosted
ORDER BY trial_count DESC
LIMIT 20
```

---

## **🎯 GRAPH PATTERN ANALYSIS**

### **13. Complete Graph Visualization (Small Sample)**
```cypher
MATCH path = (org:Organization)-[:SPONSORS]->(t:Trial)-[:CONDUCTED_AT]->(site:Site)
RETURN path
LIMIT 5
```
**Note:** In Neo4j Browser, this will show a visual graph of relationships.

### **14. Network Centrality - Most Connected Organizations**
```cypher
MATCH (o:Organization)-[r:SPONSORS|COLLABORATES_ON]->(t:Trial)
WITH o, count(r) AS connections
RETURN 
  o.name AS Organization,
  connections AS Total_Connections
ORDER BY connections DESC
LIMIT 10
```

### **15. Find Trials with Most Sites**
```cypher
MATCH (t:Trial)-[:CONDUCTED_AT]->(s:Site)
WITH t, count(s) AS site_count
RETURN 
  t.nct_id AS NCT_ID,
  t.title AS Title,
  site_count AS Site_Count
ORDER BY site_count DESC
LIMIT 10
```

---

## **🔍 SEARCH & FILTER QUERIES**

### **16. Search Trials by Title (Keyword)**
```cypher
MATCH (t:Trial)
WHERE t.title CONTAINS 'ovarian' OR t.title CONTAINS 'cancer'
RETURN 
  t.nct_id AS NCT_ID,
  t.title AS Title,
  t.status AS Status
LIMIT 20
```

### **17. Find Trials by Status**
```cypher
MATCH (t:Trial {status: 'RECRUITING'})
RETURN 
  t.nct_id AS NCT_ID,
  t.title AS Title,
  t.phase AS Phase
LIMIT 20
```

### **18. Find Organizations by Type**
```cypher
MATCH (o:Organization)
WHERE o.type = 'Sponsor'
RETURN 
  o.name AS Organization_Name,
  o.type AS Type
LIMIT 20
```

---

## **📈 AGGREGATION & STATISTICS**

### **19. Trial Distribution by Phase and Status**
```cypher
MATCH (t:Trial)
RETURN 
  t.phase AS Phase,
  t.status AS Status,
  count(t) AS Count
ORDER BY Phase, Status
```

### **20. Sites per State (Geographic Distribution)**
```cypher
MATCH (s:Site)
WHERE s.state IS NOT NULL
RETURN 
  s.state AS State,
  count(s) AS Site_Count
ORDER BY Site_Count DESC
```

### **21. Organization Type Distribution**
```cypher
MATCH (o:Organization)
RETURN 
  o.type AS Organization_Type,
  count(o) AS Count
ORDER BY Count DESC
```

---

## **🔗 COMPLEX RELATIONSHIP QUERIES**

### **22. Find All Paths Between Organization and Sites**
```cypher
MATCH path = (org:Organization)-[:SPONSORS]->(t:Trial)-[:CONDUCTED_AT]->(site:Site)
WHERE org.name CONTAINS 'National'
RETURN 
  org.name AS Sponsor,
  t.nct_id AS Trial,
  site.facility AS Site,
  length(path) AS Path_Length
LIMIT 10
```

### **23. Trials with Both Sponsor and Collaborator**
```cypher
MATCH (sponsor:Organization)-[:SPONSORS]->(t:Trial)
MATCH (collab:Organization)-[:COLLABORATES_ON]->(t)
RETURN 
  t.nct_id AS NCT_ID,
  t.title AS Title,
  sponsor.name AS Sponsor,
  collab.name AS Collaborator
LIMIT 10
```

---

## **🧪 QUICK DIAGNOSTIC QUERIES**

### **24. Check Database Health**
```cypher
CALL db.labels() YIELD label
CALL apoc.cypher.run('MATCH (n:' + label + ') RETURN count(n) as count', {}) YIELD value
RETURN label, value.count AS node_count
ORDER BY node_count DESC
```

### **25. List All Relationship Types**
```cypher
CALL db.relationshipTypes() YIELD relationshipType
CALL apoc.cypher.run('MATCH ()-[r:' + relationshipType + ']->() RETURN count(r) as count', {}) YIELD value
RETURN relationshipType, value.count AS relationship_count
ORDER BY relationship_count DESC
```

### **26. Schema Overview**
```cypher
CALL db.schema.visualization()
```

---

## **🎯 USEFUL PATTERNS FOR TESTING**

### **27. Find Trials with Missing Relationships**
```cypher
MATCH (t:Trial)
WHERE NOT EXISTS {
  (t)<-[:SPONSORS]-()
} OR NOT EXISTS {
  (t)-[:CONDUCTED_AT]->()
}
RETURN 
  t.nct_id AS NCT_ID,
  t.title AS Title
LIMIT 10
```

### **28. Isolated Nodes (No Relationships)**
```cypher
MATCH (n)
WHERE NOT (n)--()
RETURN labels(n)[0] AS Node_Type, count(n) AS Count
```

### **29. Average Sites per Trial**
```cypher
MATCH (t:Trial)-[:CONDUCTED_AT]->(s:Site)
WITH t, count(s) AS site_count
RETURN 
  avg(site_count) AS Avg_Sites_Per_Trial,
  max(site_count) AS Max_Sites,
  min(site_count) AS Min_Sites
```

---

## **🚀 HOW TO RUN THESE QUERIES**

### **Option 1: Neo4j Browser**
1. Open Neo4j Browser: `https://9669e5f3.databases.neo4j.io`
2. Connect with credentials:
   - Username: `neo4j`
   - Password: `ShMc7KhfeBfeHhvZMPFbCNiie0aeUY7eMgWSg2Y7lZ8`
3. Paste query in top panel
4. Click "Run" or press `Ctrl+Enter`

### **Option 2: Neo4j Desktop / CLI**
```bash
cypher-shell -a neo4j+s://9669e5f3.databases.neo4j.io -u neo4j -p 'ShMc7KhfeBfeHhvZMPFbCNiie0aeUY7eMgWSg2Y7lZ8'
# Then paste query
```

### **Option 3: Python Script**
```python
from neo4j import GraphDatabase

driver = GraphDatabase.driver(
    "neo4j+s://9669e5f3.databases.neo4j.io",
    auth=("neo4j", "ShMc7KhfeBfeHhvZMPFbCNiie0aeUY7eMgWSg2Y7lZ8")
)

with driver.session() as session:
    result = session.run("MATCH (t:Trial) RETURN count(t) AS count")
    print(result.single()["count"])
```

---

## **📝 EXPECTED RESULTS**

Based on current database:
- **Trials:** ~30
- **Organizations:** ~37
- **Sites:** ~860
- **PIs:** 0 (extraction needs enhancement)
- **Relationships:** 
  - SPONSORS: ~42
  - CONDUCTED_AT: ~868
  - LEADS: 0 (pending PI nodes)

---

**💡 Tip:** Start with queries #1 and #2 for overview, then explore specific patterns with the others!









