import sys
from .transit_planner import main
from .menu import menu

if __name__ == "__main__":
    if len(sys.argv) > 1:
        main()
    else:
        menu()
