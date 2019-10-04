import glob
import logging
import os
import uuid
import subprocess

import numpy
import pandas
from pkg_resources import resource_filename

logger = logging.getLogger(__name__)
BITEPI_BINARY = 'BitEpi.o'
OUTPUT_SUFFIXES = {
    'Purity.0.csv': 'p1',
    'Purity.1.csv': 'p2',
    'Purity.2.csv': 'p3',
    'Purity.3.csv': 'p4',
    'IG.0.csv': 'ig1',
    'IG.1.csv': 'ig2',
    'IG.2.csv': 'ig3',
    'IG.3.csv': 'ig4',
    'bestIG.csv': 'best_ig',
}


class ReturnCodeError(Exception):
    """Exception raised when BitEpi process returns non-zero."""
    pass


class Epistasis(object):
    """Calculate epistasis interactions given cases and controls."""
    def __init__(self, genotype_array, sample_array, working_directory='/tmp',
                 strict_intersect=False):
        """ Create a record of the cases and genotypes of the samples.

        :param genotype_array: An array using sample names for the
            column names, with the first column being the SNPs. The name
            for first column is not used. Each row consists of the SNP
            name, then either a 0, 1 or a 2 in each column to designate
            whether that sample is homozygous reference, heterozygous or
            homozygous alternate respectively at that SNP. Can be a list
            of lists, a numpy array (in which case the sample names will
            be the dtype names), or a pandas dataframe.
        :param sample_array: An array with the first column being sample
            names and the second column being a 1 for case, or 0 for
            control. Column names, if present, are ignored. Can be a
            list of lists, numpy array, or a pandas dataframe.
        :param working_directory: The directory into which the CSVs used
            for communication with the binary file will be stored. These
            are not deleted after use.
        :param strict_intersect: A boolean to control whether a
            ValueError should be raised if the samples in genotype_array
            don't exactly match the samples in sample_array. If False,
            analysis will be done on the intersect of the two arrays.
        :raises ValueError: If either sample_array or genotype_array are
            malformed, or don't contain matching samples and
            strict_intersect is True.
        """
        self._working_directory = working_directory
        self._genotype_list = self._convert_to_list(genotype_array)
        logger.debug("Converted genotype array:\n%s\nInto:\n%s",
                     genotype_array, self._genotype_list)
        self._sample_list = self._convert_to_list(sample_array, headers=False)
        logger.debug("Converted sample array:\n%s\nInto:\n%s", sample_array,
                     self._sample_list)
        self._validate_arrays(strict_intersect)
        self._array_csv = self._get_random_filename()
        self._array_list = self._create_array_list()
        self._write_to_csv()

    def compute_epistasis(self, p1=None, p2=None, p3=None, p4=None, ig1=None,
                          ig2=None, ig3=None, ig4=None, threads=2, sort=False,
                          best_ig=False):
        """Compute the epistasis interactions for each SNP combination.

        Call the BitEpi binary object with the provided arguments and
        return the result CSVs in a dictionary.

        Contains logging. For more verbose output set logging to
        logging.INFO or logging.DEBUG

        :param p1: Float that specifies that gini purity should be
            calculated  when cases and controls are split over
            individual SNPs. 0 <= p1 < 1, and represents the threshold
            purity that must be reached in order to be recorded. If
            p1 == 0 all SNPs will be recorded. If p == -1, all SNPs
            will be calculated, but none will be recorded (for
            benchmarking purposes). Will produce an output under
            "p1".
        :param p2: Same as p1, except pairs of SNPs are calculated
            instead of individual SNPs. Note that this is O(n^2). Will
            produce an output under "p2".
        :param p3: Same as p1, except triplets of SNPs are calculated
            instead of individual SNPs. Note that this is O(n^3). Will
            produce an output under "p3".
        :param p4: Same as p1, except quadlets of SNPs are calculated
            instead of individual SNPs. Note that this is O(n^4), and
            can take a long time, especially if p4 == 0. Will produce an
            output under "p4"
        :param ig1: Float that specifies that information gain should
            be calculated  when cases and controls are split over
            individual SNPs. In this case information gain is defined as
            ign = pn - max(pn-1 for all SNP subsets of length n), where
            p0 is defined of the purity of the population. 0 <= ig1 < 1,
            and represents the threshold information gain that must be
            reached in order to be recorded. If ig1 == 0 all SNPs will
            be recorded. if ig1 == -1, all SNPs will be calculated, but
            none will be recorded (for benchmarking purposes). Will
            produce an output under "ig1".
        :param ig2: Same as ig1, except pairs of SNPs are calculated
            instead of individual SNPs. Note that this is O(n^2). Will
            produce an output under "ig2".
        :param ig3: Same as ig1, except triplets of SNPs are calculated
            instead of individual SNPs. Note that this  isO(n^3). Will
            produce an output under "ig3".
        :param ig4: Same as ig1, except quadlets of SNPs are calculated
            instead of individual SNPs. Note that this is O(n^4), and
            can take a long time, especially if ig4 == 0. Will produce
            an output under "ig4".
        :param threads: The number of threads to use in the binary
        :param sort: Whether each result should be sorted from highest
            to lowest purity.
        :param best_ig: Creates a single output "best_ig", that contains
            only the most informative pair, triplet and quadlet for each
            SNP. This ignores the threshold arguments and takes O(n^4)
            time.
        :return:
            A dictionary of pandas dataframes, one for each output.
            Each row of a dataframe represents an interaction, except
            "best_ig", in which they represent a SNP, with
            columns:
                SNP/PAIR/TRIPLET/QUADLET_P - Gini purity using that
                    combination. From "p1/2/3/4", and all present
                    in "best_ig".
                SNP/PAIR/TRIPLET/QUADLET_IG - Information gain using
                    that combination. From "ig1/2/3/4", and all present
                    in "best_ig"
                SNP_A/B/C/D - SNPs participating in the interaction.
                    Number present depends on order of p/ig. Present in
                    "p/ig1/2/3/4
                SNP - Name of the snp, only present in "best_ig".
                PAIR/TRIPLET/QUADLET_1/2/3 - Additional SNPs
                    participating in the relevant interactions with SNP.
                    Only present in "best_ig".
        :raises bitepi.ReturnCodeError: If the binary returns a non-zero
            error code.
        :raises subprocess.TimeoutExpired: If the binary doesn't stop
            after output is complete.
        :raises ValueError: If the thresholds are set to values other
            than -1, or in the half-open range [0, 1). If threads is not
            a positive integer.
        """
        # Check threads argument
        if int(threads) != threads:
            logger.error("Got invalid argument threads=%s", threads)
            raise ValueError("threads must be a positive integer, got"
                             f" {threads}")

        binary = resource_filename(__name__, BITEPI_BINARY)
        args = [binary]
        thresholds = {
            '-p1': p1,
            '-p2': p2,
            '-p3': p3,
            '-p4': p4,
            '-ig1': ig1,
            '-ig2': ig2,
            '-ig3': ig3,
            '-ig4': ig4,
        }
        for threshold_name, value in thresholds.items():
            if value is not None:
                if value == -1:  # Benchmark only - no results
                    args.append(threshold_name)
                elif 0 <= value < 1:
                    args += [threshold_name, str(value)]
                else:
                    logger.error("Got invalid argument %s=%s", threshold_name,
                                 value)
                    raise ValueError("Thresholds (p/ig1-4) must be in the"
                                     " range [0, 1), or -1 for benchmarking."
                                     f"{threshold_name[1:]} is {value}.")
        output_prefix = os.path.join(self._working_directory, uuid.uuid4().hex)
        args += [
            '-i', self._array_csv,
            '-o', output_prefix,
            '-t', str(threads),
        ]
        if sort:
            args.append('-sort')
        if best_ig:
            args.append('-bestIG')
        logger.info("Calling: %s", ' '.join(args))
        process = subprocess.Popen(
            args=args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        with process.stdout:
            for line in process.stdout:
                logger.debug(line.decode().rstrip('\n'))
        try:
            return_code = process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            logger.error("Error when calling binary, binary is not"
                         " responding.")
            raise
        if return_code != 0:
            logger.error("Error when calling binary, got return-code %s.",
                         return_code)
            raise ReturnCodeError(
                f"BitEpi.o returned non-zero error code {return_code}")
        response_dict = {}
        output_prefix_length = len(output_prefix) + 1
        for file_name in glob.glob(f'{output_prefix}*'):
            file_suffix = file_name[output_prefix_length:]
            key = OUTPUT_SUFFIXES[file_suffix]
            logger.debug("Ingesting %s into '%s'", file_name, key)
            file_size = os.path.getsize(file_name)
            logger.debug("%s size is %s bytes", file_name, file_size)
            if file_size > 10 ** 9:
                logger.warning("%s is %s GB, processing could take a while.",
                               file_size / 10**9)
            response_dict[key] = pandas.read_csv(file_name)
        return response_dict

    def _create_array_list(self):
        """Create a single array out of genotypes and samples."""
        sample_case = dict(self._sample_list)
        preserved_indexes = [0]  # Always include the first column
        combined_samples = ['']  # Don't care about first column name
        for i, sample in enumerate(self._genotype_list[0][1:]):
            case = sample_case.get(sample)
            if case is not None:
                combined_samples.append(case)
                preserved_indexes.append(i + 1)
        array_list = [combined_samples] + [
            [genotype[i] for i in preserved_indexes]
            for genotype in self._genotype_list[1:]
        ]
        return array_list

    def _get_random_filename(self):
        """Construct a random filename in the working directory."""
        return os.path.join(self._working_directory, f'{uuid.uuid4().hex}.csv')

    def _validate_arrays(self, strict_intersect=False):
        """Ensure the array has the correct format and values."""
        sample_list = self._sample_list
        sample_list_len = len(sample_list)
        logger.debug("Validating arrays")
        if sample_list_len == 0:
            raise ValueError("sample array must contain at least one row.")
        if not all(isinstance(row, (list, tuple)) for row in sample_list):
            raise ValueError("sample array must be reducible to a list of"
                             " lists.")
        if not all(len(row) == 2 for row in sample_list):
            raise ValueError("sample array rows must each contain two"
                             " elements.")
        sample_samples = set(row[0] for row in sample_list)
        if not len(sample_samples) == sample_list_len:
            raise ValueError("First element of each row in sample array"
                             " must be unique.")
        if not all(row[1] in (0, 1) for row in sample_list):
            raise ValueError("Second element of each row in sample array must"
                             " be 0 or 1.")
        genotype_list = self._genotype_list
        genotype_list_len = len(genotype_list)
        if genotype_list_len == 0:
            raise ValueError("genotype array must contain at least one row.")
        if not all(isinstance(row, (list, tuple)) for row in genotype_list):
            raise ValueError("genotype array must be reducible to a list of"
                             " lists.")
        genotype_columns = len(genotype_list[0])
        if not genotype_columns > 1:
            raise ValueError("genotype array rows must each contain more than"
                             " one element.")
        genotype_samples = set(genotype_list[0][1:])
        if not len(genotype_samples) == genotype_columns - 1:
            raise ValueError("Each element in the first row of genotype array"
                             " must be unique.")
        if not len(set(row[0] for row in genotype_list)) == len(genotype_list):
            raise ValueError("First element of each row in genotype array"
                             " must be unique.")
        if not all(all(e in (0, 1, 2) for e in row[1:])
                   for row in genotype_list[1:]):
            raise ValueError("All elements except those in the first row and"
                             " column of genotype array must be 0, 1 or 2.")
        # Check the two lists have the same samples
        num_shared_samples = len(sample_samples & genotype_samples)
        missing_genotype_samples = len(genotype_samples) - num_shared_samples
        missing_sample_samples = len(sample_samples) - num_shared_samples
        if missing_genotype_samples:
            msg = (f"{missing_genotype_samples} samples from genotype array"
                   " are missing in sample array.")
            if strict_intersect:
                raise ValueError(msg)
            else:
                logger.warning(msg)
        if missing_sample_samples:
            msg = (f"{missing_sample_samples} samples from sample array"
                   " are missing in genotype array.")
            if strict_intersect:
                raise ValueError(msg)
            else:
                logger.warning(msg)

    def _write_to_csv(self):
        """Write an array to a csv file."""
        logger.debug(f"Writing combined array to {self._array_csv}")
        with open(self._array_csv, 'w') as output_file:
            output_file.writelines((','.join(str(v) for v in line) + '\n')
                                   for line in self._array_list)

    @classmethod
    def _convert_to_list(cls, array, headers=True):
        """Convert an array to a list or tuple."""
        if isinstance(array, numpy.ndarray):
            return ([array.dtype.names] if headers else []) + array.tolist()
        elif isinstance(array, pandas.DataFrame):
            return (([array.columns.tolist()] if headers else [])
                    + array.values.tolist())
        elif isinstance(array, (list, tuple)):
            # The header is probably not included here if we don't need it
            return array
        else:
            raise ValueError("array must be a pandas DataFrame, numpy ndarray,"
                             " list or tuple.")


if __name__ == '__main__':
    # simple example on three different input formats
    logging.basicConfig(level=logging.DEBUG)
    # logger.setLevel(logging.DEBUG)
    genotypes = [
        ['SNP', 'S1', 'S2', 'S3', 'S4'],
        ['a', 1, 2, 0, 0],
        ['b', 1, 1, 1, 1],
        ['c', 0, 1, 0, 2]
    ]
    samples = [
        ['S1', 0],
        ['S2', 1],
        ['S3', 0],
        ['S4', 0]
    ]
    epistasis = Epistasis(genotypes, samples)
    kwargs = {
        'ig2': 0,
        'p1': 0,
        'best_ig': True,
        'sort': True,
    }
    output = 'bestIG'
    result = epistasis.compute_epistasis(**kwargs)[output]
    logger.debug('List of Lists: \n %s', result)

    dtype = [(genotypes[0][0], '<U32')] + [(h, 'u1') for h in genotypes[0][1:]]
    logger.debug("genotypes dtype: %s", dtype)
    genotypes_np = numpy.array(
        [tuple(row) for row in genotypes[1:]],
        dtype=dtype
    )
    samples_np = numpy.array(
        [tuple(row) for row in samples],
        dtype=[('sample name', '<U32'), ('case', 'u1')]
    )
    epistasis_np = Epistasis(genotypes_np, samples_np)
    result_np = epistasis_np.compute_epistasis(**kwargs)[output]
    logger.debug('Numpy ndarray: \n %s', result_np)
    assert result.equals(result_np)

    genotypes_pd = pandas.DataFrame(
        genotypes[1:],
        columns=genotypes[0]
    )
    samples_pd = pandas.DataFrame(samples)
    epistasis_pd = Epistasis(genotypes_pd, samples_pd)
    result_pd = epistasis_pd.compute_epistasis(**kwargs)[output]
    logger.debug('pandas dataframe: \n %s', result_pd)
    assert result.equals(result_pd)
