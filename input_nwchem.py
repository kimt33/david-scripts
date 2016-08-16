from input_generator import InputGenerator

class NwchemInput(InputGenerator):
    def qualitycheck_keywords(self):
        """ Checks that the given keywords are appropriate
        """
        if 'units' not in self.keywords:
            self.keywords['units'] = 'angstroms'
        assert not ('hf' in self.keywords and 'dft' in self.keywords), 'You must select one of hf or dft calculation')

    def write_input(self, filename='')
        """ Writes the input in the appropriate format
        """
        if filename == '':
            filename = self.filename
        with open(filename, 'w') as f:
            f.write('start {0}\n'.format(self.basename))
            f.write('\n')
            # charge
            f.write('charge {0}\n'.format(self.charge))
            f.write('\n')
            # geometry
            f.write('geometry units {0}\n'.format(self.keywords['units']))
            for atom, coords in zip(self.atoms, self.coords):
                f.write('  {0:<3}   {1:<5.2f}f}{2:<5.2f}f}{3:<5.2f}f}\n'
                        ''.format(atom, *coords))
            f.write('end\n')
            f.write('\n')
            # basis
            f.write('basis\n')
            for atom in self.atoms:
                f.write('  {0:<3} library {1}\n'.format(atom, self.basis))
            f.write('end\n')
            f.write('\n')
            # scf
            scf_method = ''
            #hf
            if 'hf' in self.keywords:
                f.write('scf\n')
                f.write('  {0}; NOPEN {1}\n'.format(self.keywords['hf'], self.spin)
                f.write('  print low\n')
                f.write('end\n')
                f.write('\n')
                f.write('title "{0}"'.format(self.title, self.
            # dft
            if 'dft' in self.keywords:
                raise NotImplementedError
                f.write('dft\n')
                #f.write('  {0}; NOPEN {1}\n'.format(self.keywords['hf'], self.spin)
                f.write('  print low\n')
                f.write('end\n')

            # mp2

