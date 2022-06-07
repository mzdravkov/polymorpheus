"""create tasks

Revision ID: 4aa2e00d7805
Revises: 
Create Date: 2022-06-07 19:26:42.576326

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '4aa2e00d7805'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.create_table('tasks',
    sa.Column('id', sa.String(length=40), nullable=False),
    sa.Column('file', sa.String(length=1000), nullable=False),
    sa.Column('created_at', sa.String(length=24), nullable=True),
    sa.PrimaryKeyConstraint('id')
    )


def downgrade():
    op.drop_table('tasks')
