from pymongo import MongoClient

import config


class Database:
    _client = None  # shared client instance
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            # print('creating new db connection instance')
            cls._instance = super().__new__(cls)
            cls._instance._client = cls.get_client()
        return cls._instance

    @classmethod
    def get_client(cls):
        if cls._client is None:
            cls._client = MongoClient(
                config.MONGO_CONNECTION_URL,
                maxPoolSize=10,
                waitQueueMultiple=4
            )
        return cls._client

    def __init__(self):
        self._client = Database.get_client()
        self._db = self._client[config.MONGO_DB_NAME]
        self.id = self._client.options
        self._assemblies = self._db['assemblies']
        self._changes = self._db['changes']
        self._counters = self._db['counters']
        self._features = self._db['features']
        self._files = self._db['files']
        self._refseqchunks = self._db['refseqchunks']
        self._refseqs = self._db['refseqs']

    @property
    def assemblies(self):
        return self._assemblies

    @property
    def changes(self):
        return self._changes

    @property
    def counters(self):
        return self._counters

    @property
    def features(self):
        return self._features

    @property
    def files(self):
        return self._files

    @property
    def refseqchunks(self):
        return self._refseqchunks

    @property
    def refseqs(self):
        return self._refseqs
