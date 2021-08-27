import os
import argparse
from classes import *
# from likelihood.calculate_likelihood import optimize

if os.getcwd().startswith('/Users/davidmurphy'):
    ii = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/result/init_files/' \
         'YRI.nonhuman.BS1.6.CS0.0.initial.txt'
    w = '/Users/davidmurphy/GoogleDrive/linked_selection/data/pyLS/cons/wigz/nonhuman/chr1.nonhuman.wig.gz'
    # argv = 'caller.py datastruct compress_data -c chr17 -k bs_annos=(\'nonhuman.top5.cons\',)'.split()
    # argv = 'caller.py datastruct makeraw'.split()
    argv = 'caller.py mapstruct initial_stats -i {} -g -v -u'.format(ii).split()
    # argv = 'caller.py mapstruct optimize -i {} -u'.format(ii).split()

else:
    from sys import argv


__author__ = 'davidmurphy'


def main():
    program_choices = ['datastruct', 'mapstruct']
    choices = ' '.join('<{}>'.format(p) for p in program_choices)
    if len(argv) < 2:
        print 'error: caller: no program was specified'
        print 'usage: caller {}'.format(choices)
        exit(1)

    # store universal arguments in parent
    parent = argparse.ArgumentParser('caller', add_help=False)
    parent.add_argument('-i', '--init',
                        help='RunStruct initialization filepath')
    parent.add_argument('-k', '--kwargs',
                        help='"|"-delimited RunStruct kwargs, e.g. '
                             'cons=chimp|ncons=mammal|bs_annos=("a1", "a2")')
    parent.add_argument('-d', '--show_documentation', action='store_true',
                        default=False,
                        help='show documentation for a specified function')
    parent.add_argument('-c', '--chrom',
                        default=None,
                        help='specify the chromosome to be used')

    if argv[1] == 'datastruct':
        cls = datastruct.DataStruct
        function_choices = ['substitution_mask', 'conservation_mask',
                            'neutrality_mask', 'neutral_snps', 'snpcount',
                            'makeraw', 'calc_bkgd', 'compress_predictions',
                            'compress_data']
        parser = argparse.ArgumentParser(prog='datastruct',
                                         description=cls.__doc__,
                                         parents=[parent])
        parser.add_argument('class_function',
                            help='specify a function call to DataStruct',
                            choices=function_choices)
        parser.add_argument('-a', '--axt_file',
                            help='specify non-default human:outgroup aligment')
        parser.add_argument('-w', '--fwig',
                            help='specify non-default wig_file')
        parser.add_argument('-m', '--msk_file',
                            help='specify non-default neutmask file')
        parser.add_argument('-p', '--pct',
                            help='specify non-default percent cons')
        parser.add_argument('-t', '--coef',
                            help='specify selection coefficient for calc_bkgd')
        parser.add_argument('-n', '--bsanno',
                            help='specify bsanno for calc_bkgd')
        parser.add_argument('-s', '--save_file', help='save function result',
                            action='store_true')
        parser.add_argument('-x', '--del_cfg', help='delete calc_bkgd config',
                            action='store_true')

    else:
        assert (argv[1] == 'mapstruct')
        cls = mapstruct.MapStruct
        function_choices = ['optimize', 'initial_stats', 'final_stats']
        parser = argparse.ArgumentParser(parents=[parent])
        parser.add_argument('class_function',
                            help='specify a function call to MapStruct',
                            choices=function_choices)
        parser.add_argument('-o', '--pool', default=False,
                            help='multi-threading option flag',
                            action='store_true')

    # parse remaining system args
    args = parser.parse_args(argv[2:])

    # if "show_documentation" is used, print the function docs and exit
    if args.show_documentation:
        print getattr(cls, args.class_function).__doc__
        exit(1)

    # if "makeraw" command is used, execute and exit without further parsing
    if args.class_function == 'makeraw':
        runstruct.RunStruct(init=args.init, make=(args.init is None)).makeraw()
        exit(1)

    # initialize cls types and delete consumed arguments
    if cls.__name__ == 'MapStruct':
        # MapStruct-derived (several extra keywords)
        rst = cls(init=args.init, make=(args.init is None), pool=args.pool)
        del args.pool
    else:
        # ChromStruct-derived (no extra keywords)
        rst = cls(init=args.init, make=(args.init is None), chrom=args.chrom)
        # if chrom is not set, exit with error
        if args.chrom is None:
            parser.print_usage()
            message = '{}: error: chrom not specified\n'.format(parser.prog)
            parser.exit(message=message)
    del (args.init, args.chrom)

    # update attributes from kwargs and reset cls (if applicable)
    if args.kwargs:
        for kv in args.kwargs.split('|'):
            k, v = kv.split('=')
            rst[k] = eval(v)
        rst.reset()
    del args.kwargs

    # get the class-instance-bound function
    func = getattr(rst, args.class_function)
    del args.class_function

    # convert non-null arguments into keyword dict and pass to function
    kwargs = dict((k, v) for (k, v) in vars(args).items() if v)
    func(**kwargs)


if __name__ == '__main__':
    main()
