# Any Python script in the scripts folder will be able to import from this module.
import os
import os.path
import sys
import re
import traceback
from contextlib import contextmanager

LOG = None

# @contextmanager
# def snakemake_log(snakemake):
# 	global LOG
# 	try:
# 		if len(snakemake.log) > 0:
# 			LOG = open(snakemake.log[0], 'w')
# 			sys.stdout = LOG
# 			sys.stderr = LOG
# 			yield LOG
# 		else: yield sys.stdout
# 		# yield sys.stdout
# 	except Exception as err:
# 		# print("Shoot!")
# 		# traceback.print_exception(etype=None, value=err, tb=True)
# 		traceback.print_exc(file=LOG)
# 		if LOG is not None:
# 			LOG.flush()
# 			os.fsync(LOG)
# 		raise err
# 	finally:
# 		if LOG is not None:
# 			LOG.flush()
# 			os.fsync(LOG)
# 			LOG.close()

def snakemake_log_end():
	global LOG
	if LOG is not None:
		LOG.close()
		LOG = None


@contextmanager
def snakemake_log(snakemake):
	global LOG
	try:
		if len(snakemake.log) > 0:
			LOG = open(snakemake.log[0], 'w')
			orig_stdout = sys.stdout.fileno()
			orig_stderr = sys.stderr.fileno()

			dup_stdout = os.dup(orig_stdout)
			dup_stderr = os.dup(orig_stderr)

			os.dup2(LOG.fileno(), orig_stdout, inheritable=True)
			os.dup2(LOG.fileno(), orig_stderr, inheritable=True)

			sys.stdout = LOG
			sys.stderr = LOG
			
			yield LOG

		else: yield sys.stdout
		# yield sys.stdout
	except Exception as err:
		# print("Shoot!")
		# traceback.print_exception(etype=None, value=err, tb=True)
		traceback.print_exc(file=LOG)
		if LOG is not None:
			LOG.flush()
			os.fsync(LOG)
		raise err
	finally:
		if LOG is not None:
			os.close(orig_stdout)
			os.close(orig_stderr)

			os.dup2(dup_stdout, orig_stdout, inheritable=True)
			os.dup2(dup_stderr, orig_stderr, inheritable=True)

			os.close(dup_stdout)
			os.close(dup_stderr)

			LOG.flush()
			os.fsync(LOG)
			LOG.close()


def get_snakemake_param(snakemake, key, default):
	if key in snakemake.params.keys():
		return snakemake.params[key]
	else:
		return default


def guids_to_samples(SAMPLES, RUNS, sample_id_col='ID', guid_col='guid', run_id_col='run_id'):
	import pandas as pd

	new_metadatas = []
	for run_id in RUNS.ID:
		metadata = SAMPLES[[sample_id_col]].copy()
		metadata[run_id_col] = str(run_id)
		metadata[guid_col] = metadata[sample_id_col].apply(make_guid, args=(run_id,))
		new_metadatas.append(metadata)
	return pd.concat(new_metadatas)

def expand_metadata_by_runs(SAMPLES, RUNS):
	import pandas as pd

	guid_to_sample_map = guids_to_samples(SAMPLES, RUNS)
	return pd.merge(guid_to_sample_map, SAMPLES, on='ID', how='left')

	# new_metadatas = []
	# for run_id in RUNS.ID:
	# 	metadata = SAMPLES.copy()
	# 	metadata['run_id'] = str(run_id)
	# 	metadata['guid'] = metadata.ID.apply(make_guid, args=(run_id,))
	# 	new_metadatas.append(metadata)
	# return pd.concat(new_metadatas)

def make_guid(sample_id, run_id):
	return f'G{run_id}_{sample_id}'

guid_pattern = re.compile(r"G(?P<run>[\w\-]+)_(?P<sample>[\w\-]+)")
def guid_to_sample_run(guid):
	m = guid_pattern.match(guid)
	if m is None: raise Exception(f"GUID does not match expected format: '{guid}'")
	return m.groupdict()

def read_delim_auto(f, **kwargs):
	import pandas as pd

	_, ext = os.path.splitext(f)
	ext = ext.lower()
	if ext == '.tsv':
		kwa = {"sep":"\t", **kwargs}
	else:
		kwa = {**kwargs}
	return pd.read_csv(f, **kwa)

def listify(x):
	"""if x is a single string, wraps x in a list"""
	if isinstance(x,str):
		return [x]
	else: return x

def list_to_regex(lst):
	import re
	return "(" + '|'.join(re.escape(str(l)) for l in lst) + ")"
