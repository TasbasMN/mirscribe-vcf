import sqlite3

def drop_table(db_path, table_name):
    """
    Drop a table from a SQLite database.

    Args:
        db_path (str): The path to the SQLite database file.
        table_name (str): The name of the table to drop.
    """
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # Drop the table
    c.execute(f"DROP TABLE IF EXISTS {table_name}")

    conn.commit()
    conn.close()
