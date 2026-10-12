"""Issue #78: internal legends and rigid vector crops replace outside-strip checks."""
import argparse
from test_publication_pdf_assets import verify
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('review_root');verify(parser.parse_args().review_root)
