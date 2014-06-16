from flask import *
import mysql.connector

app = Flask(__name__)
cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')

@app.route("/")
def hello():
    cursor = cnx.cursor()
    cursor.execute('SELECT sampleID, description, date FROM samples, runs where samples.runID=runs.runID')
    samples = []
    for items in cursor:
        samples.append(items)
    return render_template('show_samples.html', samples=samples)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Development Server Help')
    parser.add_argument("-d", "--debug", action="store_true", dest="debug_mode",
                  help="run in debug mode (for use with PyCharm)", default=False)
    parser.add_argument("-p", "--port", dest="port",
                  help="port of server (default:%(default)s)", type=int, default=5000)

    cmd_args = parser.parse_args()
    app_options = {"port": cmd_args.port }

    # extra_dirs = ['./static/src',]
    # extra_files = extra_dirs[:]
    # for extra_dir in extra_dirs:
    #     for dirname, dirs, files in os.walk(extra_dir):
    #         for filename in files:
    #             filename = os.path.join(dirname, filename)
    #             if os.path.isfile(filename):
    #                 extra_files.append(filename)

    if cmd_args.debug_mode:
        app_options["debug"] = True
        app_options["use_debugger"] = False
        app_options["use_reloader"] = False
        app.debug=True
    else:
        print "Running with autoreload"
        app.debug=True
        app_options["debug"] = True
        app_options["use_debugger"] = True
        app_options["use_reloader"] = True

    app.run()
