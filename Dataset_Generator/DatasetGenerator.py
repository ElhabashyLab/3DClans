import subprocess
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import os
from io import StringIO


class DatasetGenerator:
    """
    This class generates fasta files of n sequences. The sequences are grouped into n_cl clusters.
    Each cluster is defined by a seed sequence and contains very similar sequences.
    The remaining sequences of a cluster are obtained by blasting the seed sequence.
    """    
    def __init__(self, output_dir="example_files/generated_fasta"):
        self.Entrez_email = "aron.wichtner@tuebingen.mpg.de"
        self.psiblast_results_path = "Dataset_Generator/psiblast_results.xml"
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)


    def generate(self, n: int, n_cl: int, seeds: list[str], file_name) -> str:
        """
        Generates a fasta file with n sequences grouped into n_cl clusters.
        :param n: number of sequences
        :param n_cl: number of clusters
        :param seeds: list of seed sequences (uids)
        :param file_name: name of the output fasta file
        :return: path to the generated fasta file
        """
        self._check_for_correct_input(n, n_cl, seeds)
        fasta_records = []
        for i in range(n_cl):
            seed_uid = seeds[i]
            seed_record = self._download_seed_record(seed_uid, i)
            seed_record_cleaned = self._clean_fasta_header(seed_record)
            fasta_records.append(seed_record_cleaned)
            # blast for n/n_cl - 1 similar sequences
            print(f"Finding similar sequences for cluster {i + 1}...")
            records_of_similar_sequences = self.get_records_of_similar_sequences(seed_record_cleaned, n // n_cl - 1, i)
            fasta_records.extend(records_of_similar_sequences)
        # save the generated fasta content to a file in self.output_dir
        output_path = f"{self.output_dir}/{file_name}"
        print(f"Saving fasta file at {output_path}")
        SeqIO.write(fasta_records, output_path, "fasta")
        return output_path
    
    
    def _download_seed_record(self, seed_uid, i):
        """
        Downloads the seed sequence from NCBI.
        :param seed_uid: UID of the seed sequence
        :param i: index of the cluster of the seed
        :return: list of SeqRecord objects
        """
        sequences = self._download_sequences([seed_uid])
        seed_record  = sequences[0]
        seed_record.id = seed_uid
        seed_record.description = f"Cluster_{i+1}_Seed"
        return seed_record
        
    
    def _check_for_correct_input(self, n, n_cl, seeds):
        """
        Checks if the input parameters are valid.
        :param n: number of sequences
        :param n_cl: number of clusters
        :param seeds: list of seed sequences (uids)
        """
        if n <= 0:
            raise ValueError("Number of sequences must be positive.")
        if n_cl > n:
            raise ValueError("Number of clusters cannot be greater than number of sequences.")
        if len(seeds) != n_cl:
            raise ValueError("Number of seed sequences must match number of clusters.")
        
    
    def _download_sequences(self, uids: list):
        """
        Downloads protein sequences from NCBI given their UIDs.
        :param uid: the UID of the sequence
        :return: the sequence as a string
        """
        Entrez.email = self.Entrez_email
        uids_as_string = ",".join(uids)
        handle = Entrez.efetch(db="protein", id=uids_as_string, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta-pearson"))
        handle.close()
        return records
    
    
    def _blast_sequence(self, record):
        """
        PSI-Blasts a protein record.
        :param record: the protein record to blast
        :return: result of PSI-BLAST
        """
        # Run PSI-BLAST with the fetched sequence
        result_handle = NCBIWWW.qblast(
            program="blastp",
            database="swissprot", # makes sure sequences have structures in alphafold DB
            sequence=str(record.seq),
            service="psi", # for normal blast use: "plain"
            word_size=3
        )
        return result_handle
    
    
    def get_records_of_similar_sequences(self, record, num_sequences, i):
        """
        Retrieves records of sequences similar to the given record sequence and makes sure the querry record is not included.
        :param record: the record to compare
        :param num_sequences: number of similar records to retrieve
        :param i: index of the cluster
        :return: list of SeqRecord objects
        """
        blast_results = self._blast_sequence(record)
        similar_records = self._parse_hits(record, num_sequences, blast_results)
        for record in similar_records:
            record.description = f"Cluster_{i+1}_Similar"
        print(f"Downloaded {len(similar_records)} of {num_sequences} similar sequences for cluster {i + 1}.")
        return similar_records


    def _parse_hits(self, seed_record, max_hits, blast_results):
        """
        Parses the PSI-BLAST hits from the XML result file and downloads the corresponding sequences from NCBI.
        It leaves out the seed sequence if it is among the hits.
        :param seed_record: the blast query record
        :param max_hits: maximum number of similar records to retrieve
        :param blast_results: the result of PSI-BLAST as a file-like object
        :return: list of SeqRecord objects not longer than max_hits
        """
        fasta_records = []
        hit_ids = []
        blast_records = list(NCBIXML.parse(blast_results))
        # get results of the last PSI iteration
        final_blast_record = blast_records[-1]
        for alignment in final_blast_record.alignments:
            accession = alignment.accession
            if accession not in hit_ids and len(hit_ids) < max_hits and accession != seed_record.id:
                hit_ids.append(accession)
            else:
                continue
        fasta_records = self._download_sequences(hit_ids)
        fasta_records_cleaned = [self._clean_fasta_header(record) for record in fasta_records]
        return fasta_records_cleaned


    def _save_to_file(self, filename, content):
        """
        Saves a file with a given content.
        :param filename: path to the file
        :param content: content to write into the file
        """
        with open(filename, "w") as f:
            f.write(content)
            
    
    def _clean_fasta_header(self, record):
        """
        Cleans the fasta header of a SeqRecord object so it only contains the UniProt ID.
        :param record: the SeqRecord object
        :return: the cleaned SeqRecord object
        """
        if "|" in record.id:
            record.id = record.id.split("|")[1].split(".")[0]
            return record
        else:
            record.id = record.id.split(".")[0]
            return record

            
# test
example_seeds = ["P68871", "Q99895", "P42212", "P00734", "P69905", "P0A6F5", "Q8N3C0", "P00519", "P00846", "P00390", "P02754", "Q8RWR1"]
generator = DatasetGenerator()
#generator.generate(15, 2, ["P42212", "Q99895"], "test.fasta")
#for i in range(4):
    #generator.generate(15, 3, example_seeds[i*3:(i+1)*3], f"example_{i+1}.fasta")
