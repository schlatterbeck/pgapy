VERSION="RELEASE".replace ('V_', '').replace ('_', '.')

if __name__ == "__main__" :
    print '#define VERSION "%s"' % VERSION
