import scripter
import os.path

class AlignmentsFilenameParser(scripter.FilenameParser):
    """
    for alignments produced by align.py
    """
    def __init__(self, filename, controls = {}, *args, **kwargs):
        if not os.path.splitext(filename)[1] == '.bam':
            raise InvalidFileException(filename)
        super(AlignmentsFilenameParser, self).__init__(filename, *args, **kwargs)
        
        sample = self.protoname
        # check controls
        if controls.has_key(sample):
            run_name, control = controls[sample]
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
            run_name = sample
            controls[sample] = (sample, None)
            
        self.run_name = run_name  
        self.output_dir = os.path.join(self.output_dir, run_name)
        self.check_output_dir(self.output_dir)