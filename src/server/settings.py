"""
Django settings for server project.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.6/ref/settings/
"""

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.6/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 't9+m%qyni5%=__s8brz#tf#lv^1wy6)zj#m_2re&(_c(!_pixl'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

TEMPLATE_DEBUG = True

ALLOWED_HOSTS = ['*']

import sys
TESTING = sys.argv[1:2] == ['test']

# for Django Celery

BROKER_URL = os.environ.get('AMQP_URL', 'amqp://guest@localhost:5672//')


if TESTING:
  CELERY_ALWAYS_EAGER = True  # skip the daemon
else:
  CELERY_SEND_TASK_SENT_EVENT = True

# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django_assets',
    'edge',
)

if TESTING:
    INSTALLED_APPS += ('django_nose',)
    TEST_RUNNER = 'django_nose.NoseTestSuiteRunner'

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    #'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'server.urls'

WSGI_APPLICATION = 'server.wsgi.application'


# Database
# https://docs.djangoproject.com/en/1.6/ref/settings/#databases

DATABASES = {
    'sqlite': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
    },
    'mysql': {
        'ENGINE': 'django.db.backends.mysql',
        'OPTIONS': { "init_command": "SET storage_engine=INNODB;" },
        'HOST': os.getenv("DB_HOST", ""),
        'PORT': "",
        'NAME': os.getenv("DB_NAME", ""),
        'USER': os.getenv("DB_USER", ""),
        'PASSWORD': os.getenv('DB_PASSWORD', ""),
        'ATOMIC_REQUESTS': True,
    }
}

DEFAULT_DB = 'mysql'

DATABASES['default'] = DATABASES[DEFAULT_DB]
if TESTING:
    other_dbs = [db for db in DATABASES if db != 'default']
    for db in other_dbs:
        del DATABASES[db]

# Internationalization
# https://docs.djangoproject.com/en/1.6/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.6/howto/static-files/

STATIC_URL = '/static/'

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        }
    },
    'handlers': {
        'console': {
            'level': 'DEBUG',
            'class': 'logging.StreamHandler'
        }
    },
    'loggers': {
        #'django.db.backends': { 'level': 'DEBUG', 'handlers': ['console'], },
    },
}


# NCBI blast
NCBI_DIR = os.getenv('NCBI_DIR', BASE_DIR + '/../ncbi')
NCBI_BIN_DIR = os.getenv('NCBI_BIN_DIR', NCBI_DIR + '/bin')
NCBI_DATA_DIR = os.getenv('NCBI_DATA_DIR', NCBI_DIR + '/blastdb')

# Primer3
PRIMER3_BIN = os.getenv('PRIMER3_BIN', BASE_DIR + '/../primer3/primer3_core')
PRIMER3_CONFIG_DIR = os.getenv('PRIMER3_CONFIG_DIR', BASE_DIR + '/../primer3/primer3_config')


NOSE_ARGS = [
  '--with-coverage',
  '--cover-package=edge',
]
