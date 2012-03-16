import scripter
import os.path
from bioplus.genometools import guess_bam_genome, NoMatchFoundError

class BAMFilenameParser(scripter.FilenameParser):
    """
    for alignments produced by align.py
    """
    def __init__(self, filename, controls = {}, *args, **kwargs):
        if not os.path.splitext(filename)[1] == '.bam':
            raise scripter.InvalidFileException(filename)
        super(BAMFilenameParser, self).__init__(filename, *args, **kwargs)
        
        sample = self.protoname
        # check controls
        if controls.has_key(sample):
            sample_name, control = controls[sample]
            scripter.debug('%s has control %s', sample, control)
            if control is None:
                self.control_file = None
            else:
                self.control_file = os.path.join(self.input_dir,
                                                 control + '.bam')
        elif sample in [v[1] for v in controls.values()]:
            scripter.debug('%s is a control, aborting', sample)
            raise scripter.InvalidFileException
        else:
            scripter.debug('%s has no control indicated, continuing anyway',
                           sample)
            # not in setup.txt, make an entry in controls
            self.control_file = None
            sample_name = sample
            controls[sample] = (sample, None)
            
        self.sample_name = sample_name
        self.output_dir = os.path.join(self.output_dir, sample_name)