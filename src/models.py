from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

db = create_engine('sqlite:///db.sqlite', echo=True)


class Task(Base):
    __tablename__ = 'tasks'
    # sha1 sum of the file
    id = Column(String(40), primary_key=True)
    # file name
    file = Column(String(1000), nullable=False)
    # task creation timestamp in ISO8601 format
    created_at = Column(String(24))