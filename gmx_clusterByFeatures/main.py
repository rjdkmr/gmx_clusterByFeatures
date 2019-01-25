import sys


def main():
    options = {'cluster':'Perform clustering using features and extract clustered trajectories',
              }

    if len(sys.argv)<=1:
        show_help(options)
        sys.exit(-1)

    if sys.argv[1] not in options:
        print(' ERROR: "{0}" is not an accepted option.\n' .format(sys.argv[1]))
        show_help(options)
        sys.exit(-1)

    if sys.argv[1] == 'cluster':
        from .gmx_clusterByFeatures import cluster
        cluster(sys.argv[1:])
        
def show_help(options):
    print(' ==============================')
    print(' Usage:')
    print(' gmx_clusterByFeatures <Option>\n')
    print(' ---------------------')
    print(' Use following options:')
    print(' ---------------------\n')

    for tool in options:
        print('\t{0} : {1}\n'.format(tool, options[tool]))

    print(' ==============================')


if __name__=="__main__":
    main()
