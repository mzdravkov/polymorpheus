import duckdb

con = duckdb.connect(database='db.duckdb', read_only=False)


con.execute(
"""
CREATE TABLE IF NOT EXISTS tasks (
	id VARCHAR(40) NOT NULL, 
	file VARCHAR(1000) NOT NULL, 
	created_at VARCHAR(24), 
	PRIMARY KEY (id)
)
"""
)

def insert_