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

		else:
			yield sys.stdout
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
