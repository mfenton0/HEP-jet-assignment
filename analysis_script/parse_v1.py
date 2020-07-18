from HEAH import parse_function
import sys, os, time
def main():

    root_file = sys.argv[1]  #your root file path
    time_tag = time.strftime("%Y_%m_%d_%H_%M", time.localtime())
    FILENAME = ""
    seq = ("parsed_", time_tag, ".txt")
    PATH = os.path.join(os.getcwd(), FILENAME.join(seq))
    parsed_event = parse_function.parse_function(root_file, PATH)
    particle_numbers = parsed_event.count_particles()
    saved_path = parsed_event.save_parsed_file()
    print("The data has been parsed! There are {0} particles in this events. \n The result has been store in '{1}'.".format(particle_numbers, saved_path))


if __name__ == "__main__": 
    main()
