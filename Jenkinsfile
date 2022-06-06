pipeline {
  agent none

  options {
    buildDiscarder(logRotator(numToKeepStr:'50'))
  }

  // TODO: checkout https://pre-commit.com/
  // https://github.com/PyCQA/prospector#pre-commit
  stages {
    // formatter
    // ---------
    // https://github.com/ambv/black
    // https://github.com/google/yapf
    // https://medium.com/3yourmind/auto-formatters-for-python-8925065f9505
    //
    // linter
    // ------
    // https://github.com/PyCQA/prospector
    // https://prospector.readthedocs.io/en/latest/supported_tools.html#optional-extras
    stage('lint') {
      agent { docker { alwaysPull true; image 'python:3' } }
      steps {
        // numpy is installed first because it is a chicken-or-the-egg dependency.
        // ASAP is installed to silence 'cannot import <dependency>' warnings.
        //
        // autoformat with black to silence warnings that can be auto-fixed.
        // black excludes .git/ and .venv/ by default.
        //
        // max-line-length default is 80
        sh '''\
          python -m venv .venv
          . .venv/bin/activate
          pip install --upgrade pip
          pip install numpy
          pip install .
          pip install 'prospector[with_everything] >=1.1.7' black
          black --quiet .
          prospector \
            --max-line-length 120 \
            --output-format pylint:prospector.log \
            --zero-exit \
            --strictness veryhigh \
            --doc-warnings \
            --test-warnings \
            --member-warnings \
            --full-pep8 \
            --with-tool pyroma \
            --with-tool vulture \
            --with-tool frosted \
            --with-tool mypy
          '''.stripIndent()
        recordIssues(tools: [pyLint(pattern: 'prospector.log')])
      }
      post {
        cleanup { deleteDir() }
      }
    }
  }
}
