import os

def export_envs(request):
    data = {}
    data['FULLSTORY_ORG_ID'] = os.getenv('FULLSTORY_ORG_ID', '')
    return data
