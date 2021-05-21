from django import forms


class FragmentForm(forms.Form):
    name = forms.CharField()
    sequence = forms.CharField()
    circular = forms.BooleanField(required=False)
