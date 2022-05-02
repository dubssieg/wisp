"""
    output = None
    if args.index != None:
        make_datasets(
            input_dir=INPUT_PATH, path=OUTPUT_PATH,job_name=JOB, db_name=DATABASE, sampling=SAMPLING, kmer_size=KMER_SIZE, func=FUNC, ratio=RATIO, read_size=RS, classif_level=taxa)
        make_model(OUTPUT_PATH, JOB, classif_level, DATABASE)
    if args.test:
        output = test_model(OUTPUT_PATH, JOB, DATABASE)
    else:
        if TAXON == 'order':
            create_unk_sample(UNK_GENOME, SAMPLING, OUTPUT_PATH,
                              JOB, kmer_size=KMER_SIZE, func=FUNC, ratio=RATIO)
            output = test_unk_sample(OUTPUT_PATH, JOB, DATABASE, TAXON)
        elif TAXON == 'family':
            create_unk_sample(UNK_GENOME, SAMPLING, OUTPUT_PATH,
                              JOB, kmer_size=KMER_SIZE, func=FUNC, ratio=RATIO)
            output = test_unk_sample(OUTPUT_PATH, JOB, DATABASE, 'order')
            make_datasets_family(
                INPUT_PATH, OUTPUT_PATH, DATABASE, SAMPLING, output['Reads summation order'], kmer_size=KMER_SIZE, func=FUNC, ratio=RATIO)
            make_model_family(OUTPUT_PATH, DATABASE,
                              output['Reads summation order'])
            create_unk_sample(UNK_GENOME, SAMPLING, f"{OUTPUT_PATH}{output['Reads summation order']}/",
                              JOB, kmer_size=KMER_SIZE, func=FUNC, ratio=RATIO)
            output = {**output, **test_unk_sample(
                f"{OUTPUT_PATH}{output['Reads summation order']}/", JOB, DATABASE, 'family')}
"""

"""    
    # declaring args
    parser.add_argument(
        "genome_path", help="genome/read to evaluate", type=str)
    parser.add_argument("job_name", help="Name of job", type=str)
    parser.add_argument(
        "database_name", help="Name of database, raises error if not exists", type=str)
    parser.add_argument('classiftarget',
                        choices=['order', 'family'],
                        help='defines the target classification level', type=str)
    parser.add_argument(
        "-i", "--index", type=str, help="path to genome database")
    parser.add_argument(
        "-o", "--output", type=str, help="path to output data folder")
    parser.add_argument(
        "-f", "--func", type=str, choices=['mean'], help="fuction to filter out results")
    parser.add_argument(
        "-m", "--method", type=str, choices=['build_only', 'pred'], help="fuction to filter out results")
    parser.add_argument(
        "-t", "--test", help="tests the model with test set", action="store_true")
    parser.add_argument("-s", "--subsampling", type=int,
                        help="give a sampling number for database creation")
    parser.add_argument("-r", "--readsize", type=int,
                        help="give a sampling size for reads")
    parser.add_argument("-l", "--limiter", type=float,
                        help="factor for func")
    parser.add_argument("-k", "--kmer", type=int,
                        help="give a size for kmers at sampling ; be careful with database version")
    #parser.add_argument("-n", "--numrun", type=int,help="give a number of runs for unknown")

    # executing args
    args = parser.parse_args()
    # task_00 > name of the job (under which name they will be saved)
    JOB: str = args.job_name
    # task_00 > name of the database (under which name they will be saved)
    DATABASE: str = args.database_name
    # "/udd/sidubois/Stage/Genomes/Lactobacillales_Lactobacillus_amylovorus.fna" > full path to genome to evaluate
    UNK_GENOME: str = args.genome_path
    # -s 50 > number of reads to cut inside each genome
    SAMPLING: int = args.subsampling if args.subsampling != None else 500
    # -i "/udd/sidubois/Stage/Genomes/"  > path to genomes
    INPUT_PATH: str = args.index
    # -o "data/" > path to output, w. default
    OUTPUT_PATH: str = args.output if args.output != None else "data/"
    # -k 5 > size of indexing
    KMER_SIZE: int = args.kmer if args.kmer != None else 5
    # -f mean > regression func used to filter out results, from stats package
    FUNC: str | None = getattr(
        statistics, args.func) if args.func != None else None
    # -l 1.5 > factor applied to func to filter out results, def 1.5
    RATIO: float = args.limiter if args.limiter != None else 1.5
    # -c order > classification level to target
    TAXON: str = args.classiftarget
    # -r 10000 > size of subsampled reads
    RS: int = args.readsize if args.readsize != None else 10000

"""
