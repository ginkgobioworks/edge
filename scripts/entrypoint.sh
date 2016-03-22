#!/usr/bin/env bash
set -e

. $HOME/env/bin/activate

DJANGO_IP=${DATASTORE_IP:-0.0.0.0}
DJANGO_PORT=${DATASTORE_PORT:-8080}
MANAGE=edge_manage.py

while [[ $# -gt 0 ]]
do
    case "$1" in
        collectstatic )
            $MANAGE collectstatic --noinput
            shift 1
            ;;
        syncdb )
            $MANAGE syncdb --noinput
            shift 1
            ;;
        migrate )
            $MANAGE migrate --noinput
            shift 1
            ;;
        init )
            $MANAGE syncdb --noinput
            $MANAGE migrate --noinput
            $MANAGE collectstatic --noinput
            shift 1
            ;;
        run* )
            $MANAGE runserver $DJANGO_IP:$DJANGO_PORT
            shift 1
            ;;
        * )
            $MANAGE $1
            shift 1
            ;;
    esac
done
