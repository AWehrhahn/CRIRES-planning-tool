import datetime
from astropy.time import Time
from astropy import units as u

from .transit_planner import transit_calculation
from .interactive_graph import create_interactive_graph

class MenuController:
    def __init__(self, hello="Hello"):
        self.hello = hello

        # Parse the methods and docstrings
        methods = [m for m in dir(self) if m.startswith('do_')]
        self.methods = {m[3:]: getattr(self, m) for m in methods}

        index = sorted(list(self.methods.keys()))
        doc_strings = [f"{i}\t{self.methods[i].__doc__}"  for i in index]
        self.width = max([len(d) for d in doc_strings])
        self.menu_string = "\n".join(doc_strings)

    def do_q(self):
        """ Exit Program """
        exit()

    def generate_menu(self):
        print(self.hello)
        print("=" * self.width)
        print(self.menu_string)
        print("=" * self.width)

        while True:
            end = input("Select an option: ")
            if end not in self.methods.keys():
                print("Invalid Option %i" % end)
                continue
            else:
                break

        return self.methods[end]()

class MenuMain(MenuController):
    def __init__(self):
        super().__init__(self)
        self.hello = """
   _____ _____  _____ _____  ______  _____        _____  _                   _               _______          _ 
  / ____|  __ \|_   _|  __ \|  ____|/ ____| _    |  __ \| |                 (_)             |__   __|        | |
 | |    | |__) | | | | |__) | |__  | (___ _| |_  | |__) | | __ _ _ __  _ __  _ _ __   __ _     | | ___   ___ | |
 | |    |  _  /  | | |  _  /|  __|  \___ \_   _| |  ___/| |/ _` | '_ \| '_ \| | '_ \ / _` |    | |/ _ \ / _ \| |
 | |____| | \ \ _| |_| | \ \| |____ ____) ||_|   | |    | | (_| | | | | | | | | | | | (_| |    | | (_) | (_) | |
  \_____|_|  \_\_____|_|  \_\______|_____/       |_|    |_|\__,_|_| |_|_| |_|_|_| |_|\__, |    |_|\___/ \___/|_|
                                                                                      __/ |                     
                                                                                     |___/                      
        """
        # self.hello = "Welcome to the CRIRES+ Planning Tool"

    def do_1(self):
        """ Single Planet Calculation """
        while True:
            planet = input("Choose Planet: ")
            if planet == "":
                print("Must choose a planet")
            else:
                break
            
        date_start = input("Starting Date [today]: ")
        date_end = input("End Date [Start date + 1 year]: ")
        savefile = input("Output filename [data.csv]: ")
        plot = input("Create plot [no]: ")

        if plot == "":
            plot = False
        else:
            plot = bool(plot)
        if plot:
            plot_file = input("Plot filename [plot.html]: ")
            if plot_file == "":
                plot_file = "plot.html"

        if date_start == "":
            date_start = str(datetime.date.today())
        if date_end == "":
            date_end = Time(date_start) + 1 * u.year
            date_end = str(date_end.datetime)            
        if savefile == "":
            savefile = "data.csv"

        df = transit_calculation(planet, date_start, date_end, mode="planets")
        df.to_csv(savefile)

        if plot:
            plot_file1 = "_1.".join(plot_file.rsplit(".", 1))
            plot_file2 = "_2.".join(plot_file.rsplit(".", 1))
            create_interactive_graph(df, plot_file1, plot_file2)
        
        pass

    def do_2(self):
        """ Full Catalogue Calculation """
        criteria = input("Choose planet/star criteria [default]: ")
        date_start = input("Starting Date [today]: ")
        date_end = input("End Date [Start date + 1 year]: ")
        savefile = input("Output filename [data.csv]: ")
        plot = input("Create plot [no]: ")

        if plot == "":
            plot = False
        else:
            plot = bool(plot)
        if plot:
            plot_file = input("Plot filename [plot.html]: ")
            if plot_file == "":
                plot_file = "plot.html"


        if criteria == "":
            criteria = None
        if date_start == "":
            date_start = str(datetime.date.today())
        if date_end == "":
            date_end = Time(date_start) + 1 * u.year
            date_end = str(date_end.datetime)            
        if savefile == "":
            savefile = "data.csv"
        
        df = transit_calculation(criteria, date_start, date_end, mode="criteria")
        df.to_csv("data.csv")

        if plot:
            plot_file1 = "_1.".join(plot_file.rsplit(".", 1))
            plot_file2 = "_2.".join(plot_file.rsplit(".", 1))
            create_interactive_graph(df, plot_file1, plot_file2)

        pass

def menu():
    main_menu = MenuMain()
    main_menu.generate_menu()


if __name__ == "__main__":
    menu()
