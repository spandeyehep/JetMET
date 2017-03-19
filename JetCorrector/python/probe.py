

class Store(object):
    def __init__( self ):
        self.attributes = []
    def append_attribute( self, name):
        if name not in self.attributes:
            print "Adding attribute %r" % name
            self.attributes.append( name )

class Probe( Store ):
     
    def __init__( self ):
        super(Probe, self ).__init__()

    def __getattribute__( self, name):
        try:
            return object.__getattribute__( self, name )
        except AttributeError as e:
            super(Probe, self).append_attribute(name )
            return 1.
            raise e
            

def probing( formula ):
    probe = Probe()
    try:
        formula( probe )
    except ZeroDivisionError:
        pass

    print "found attributes %s" % ",".join(probe.attributes)

formula = lambda x: x.Px/(x.Px-x.Py)+sqrt( - x.Pz )

probing( formula )
