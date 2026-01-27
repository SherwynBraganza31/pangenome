import json, requests, os, re, time, zlib, pandas as pd
from requests.adapters import HTTPAdapter, Retry
from urllib.parse import urlparse, parse_qs, urlencode
from xml.etree import ElementTree

class keggRestHandler:
    """
            Created based on the manual provided by https://www.kegg.jp/kegg/rest/keggapi.html.
            Handles REST API interfacing and result retrieval.
            For more background information refer to the docs from www.kegg.jp

            Attributes:
            -----------
            job_id : The id of the job submitted to UniProt
            response: The decoded response in JSON format
            headers:

        """

    def __init__(self, polling_interval: int = 3, num_retries: int = 5):
        self.POLLING_INTERVAL = polling_interval  # set number of seconds to wait before retrying
        self.API_URL = "https//rest.kegg.jp/"

        self.retries = Retry(total=num_retries, backoff_factor=0.25,
                             status_forcelist=[500, 502, 503, 504])
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=self.retries))

        self.job_id = None
        self.response = None
        self.headers = None

    def check_response(self, response: requests.Response):
        """
        Checks the status of the HTTP request and raises it if there's an error.
        Static Helper function only used within the class.

          @params:
            response : requests.Response
            request.Response object containing information of the request

          @returns:
            None
        """
        try:
            response.raise_for_status()
        except requests.HTTPError:
            print(response.json())
            raise

    def generate_ec_url(self, ec_num:str):
        """
        Generates the EC encoding to be injected into the url.

        :param ec_num: The enzyme number encoding to get pathway resutls for
        :return:
        """
        if '-' in ec_num:
            raise Exception(f'Incomplete EC number given : {ec_num}')

        url_tag = f'{self.API_URL}get/ec:{ec_num}/'
        return url_tag


    def get_ec_results(self, source_dir:str):
        postproc_results_dir = source_dir + 'postprocessing_results/'

        with open(postproc_results_dir + 'ecNumbers.txt', 'r') as f:
            ec_numbers = f.read().splitlines()

        complete_ec = [x for x in ec_numbers if '-' not in x]

        for ec in complete_ec:
            response = requests.get(self.generate_ec_url(ec))
            if response.status_code != 200:
                raise Exception(f"Failed to retrieve data for {ec}; http status code - {response.status_code}")

            data = response.text

            parsed = {}
            current_key = None

            for line in data.splitlines():
                if not line.strip():
                    continue
                if line[:12].strip():  # New key
                    current_key = line[:12].strip()
                    parsed[current_key] = line[12:].strip()
                else:  # Continuation of previous key
                    parsed[current_key] += " " + line[12:].strip()

            with open(postproc_results_dir +'kegg_ec_map/' + ec + '.json', 'w') as f:
                json.dump(parsed, f)

