from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import os

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
            seed_record = self._download_seed_records(seed_uid, i)
            seed_record = self._clean_fasta_header(seed_record)
            fasta_records.append(seed_record)
            # blast for n/n_cl - 1 similar sequences
            print(f"Finding similar sequences for cluster {i + 1}...")
            records_of_similar_sequences = self.get_records_of_similar_sequences(seed_record, n // n_cl - 1, i)
            fasta_records.extend(records_of_similar_sequences)
        # save the generated fasta content to a file in self.output_dir
        output_path = f"{self.output_dir}/{file_name}"
        print(f"Saving fasta file at {output_path}")
        SeqIO.write(fasta_records, output_path, "fasta")
        return output_path
    
    
    def _download_seed_records(self, seed_uid, i):
        """
        Downloads the seed sequences from NCBI.
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
        Downloads a protein sequence from NCBI given its UID.
        :param uid: the UID of the sequence
        :return: the sequence as a string
        """
        Entrez.email = self.Entrez_email
        uids_as_string = ",".join(uids)
        handle = Entrez.efetch(db="protein", id=uids_as_string, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return records
    
    
    def _blast_sequence(self, record):
        """
        PSI-Blasts a protein record.
        :param record: the protein record to blast
        :return: result of BLAST
        """
        # Run PSI-BLAST with the fetched sequence
        result_handle = NCBIWWW.qblast(
            program="blastp",
            database="swissprot", # makes sure sequences have structures in alphafold DB
            sequence=str(record.seq),
            service="psi",
            word_size=3
        )
        return result_handle
    
    
    def get_records_of_similar_sequences(self, record, num_sequences, i):
        """
        Retrieves records of sequences similar to the given record sequence.
        :param record: the record to compare
        :param num_sequences: number of similar records to retrieve
        :param i: index of the cluster
        :return: list of SeqRecord objects
        """
        results_handle = self._blast_sequence(record)
        self._save_to_file(self.psiblast_results_path, results_handle.read())
        results_handle.close()
        similar_records = self._parse_hits(num_sequences)
        for record in similar_records:
            record.description = f"Cluster_{i+1}_Similar"
        print(f"Downloaded {len(similar_records)} of {num_sequences} similar sequences for cluster {i + 1}.")
        return similar_records


    def _parse_hits(self, max_hits):
        """
        Parses the BLAST hits from the XML result file.
        :param max_hits: maximum number of similar records to retrieve
        :return: list of SeqRecord objects not longer than max_hits
        """
        fasta_records = []
        hit_ids = []
        with open(self.psiblast_results_path) as f:
            blast_records = list(NCBIXML.parse(f))
        # get results of the last PSI iteration
        final_blast_record = blast_records[-1]
        for alignment in final_blast_record.alignments:
            accession = alignment.accession
            if accession not in hit_ids and len(hit_ids) < max_hits:
                hit_ids.append(accession)
            else:
                break
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
# other possible seeds: "Q99895", "P42212"
# seeds = ["P68871"]
# generator = DatasetGenerator()
# generator.generate(30, 1, seeds, "test.fasta")
