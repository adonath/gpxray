#!/bin/bash

# --8<--8<--8<--8<--
#
# Copyright (C) 2020 Smithsonian Astrophysical Observatory
#
# This file is part of saotrace
#
# dist is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# -->8-->8-->8-->8--

VERSION=2.0.5

trap cleanup EXIT

TMP_PFX=/tmp/run-saotrace-$$

cleanup () {
    /bin/rm -rf "$TMP_PFX.*"
}

die () {
    echo >&2 "$@"
    exit 1
}

is_in_path () {
    command -v $1 > /dev/null 2>&1 
}

LAUNCHER=${RUN_SAOTRACE_LAUNCHER:-any}
WORKDIR=$PWD
IMAGE=saotrace:${VERSION}
CUID=0
CGID=0

help () {

    cat <<EOF
usage: run-saotrace [options] <command> [options for <command>]

  Options
  -l,--launcher TEXT   specify launcher to start the container, 
                          one of 'docker', 'podman', or 'any'
  --workdir DIR        cd to DIR in the container before starting the command.
                          defaults to '${WORKDIR}'
  --image IMAGE        use the container image. defaults to '$IMAGE'.
  -h,--help            output short usage instructions and exit
  --man,--manual       output a manual page and exit
  -v,--version         output version informaiton and exit

  <command> is any command in the container, such as
    trace-nest
    plist
    pset
    bash
EOF
    exit 0
}

_roff () {
    local file=$1
    sed -n -e '
/^# <<< roff/,/^# >>> roff/{
/^# <<< roff/d
/^# >>> roff/d
p
}
' $0 > "$file"
    test $? -ne 0 && die "unable to create temporary file"
}

manual () {
    local tmpfile
    tmpfile=$(mktemp -- "$TMP_PFX.XXXXXXX")
    test $? -ne 0 && die "unable to create temporary file"
    _roff $tmpfile
    man $tmpfile
    exit 0
}

# these need to be global for the run-time generated docker shell
# functions to work.
DOCKER_EXE=
PODMAN_EXE=

setup_docker () {

    local launcher=$1

    case $launcher in

        podman )
            is_in_path podman || die "podman is not in your PATH"
            PODMAN_EXE=$(command -v podman)
            docker () { ${PODMAN_EXE} "$@" ; }
	    CUID=0
            CGID=0
            ;;
        
        docker )
            is_in_path docker || die "docker is not in your PATH"
            DOCKER_EXE=$(command -v docker)
            local in_docker_group=no
            for group in $(groups $(id --name --user))
            do
                test "$group" = "docker" && in_docker_group=yes
            done
            if [ $in_docker_group = yes ];  then
                docker () { ${DOCKER_EXE} "$@" ; }
            else
                docker () { sudo ${DOCKER_EXE} "$@" ; }
            fi

            CUID=$(id --user)
            CGID=$(id --group)       \
            ;;

        any )
            if is_in_path podman; then
                launcher=podman
            elif is_in_path docker; then
                launcher=docker
            else
                die "unable to find either docker or podman"
            fi
            setup_docker "$launcher"
            ;;

        * )
            die "Unknown container launcher: $launcher"
            ;;

    esac
}    

_docker_run () {
    local entry_point=$1
    test $# -ne 0 && shift

    docker run -i -t --rm                       \
	-v $WORKDIR:/workdir                    \
	--user ${CUID}:${CGID}                  \
        --entrypoint "$entry_point"             \
        --platform linux/amd64  \
	$IMAGE                                  \
        "$@"
}


while [ $# -ne 0 ];  do

      arg=$1

      case "$arg" in

          --launcher=* )
              LAUNCHER=$(echo "$arg" | sed -e 's/^--use//')
              ;;

          --launcher | -l )
              shift
              test $# -ne 0 || die "missing argument to $arg"
              LAUNCHER=$1
              ;;

          --workdir=* )
              WORKDIR=$(echo $1 | sed -e 's/^--workdir//')
              ;;

          --workdir )
              shift
              test $# -ne 0 || die "missing argument to $arg"
              WORKDIR=$1
              ;;

          --image=*)
              IMAGE=$(echo $1 | sed -e 's/^--image//')
              ;;

          --image)
              shift
              test $# -ne 0 || die "missing argument to $arg"
              IMAGE=$1
              ;;

          -h|--help )
              help
              ;;

          --man|--manual )
              manual
              ;;

          -v|--version)
              echo saotrace $VERSION
              exit 0
              ;;

          -x )
              set -x
              ;;

          -- )
              shift
              break
              ;;

          -* )
              die "unknown option $1"
              ;;
          *)
              break
              ;;
      esac

      shift
done


setup_docker "$LAUNCHER"
test $# -eq 0 && die "no command was specified; use the --help option for more details"
_docker_run "$@"
exit $?

# <<< roff
.\" Automatically generated by Pod::Man 4.10 (Pod::Simple 3.35)
.\"
.\" Standard preamble:
.\" ========================================================================
.de Sp \" Vertical space (when we can't use .PP)
.if t .sp .5v
.if n .sp
..
.de Vb \" Begin verbatim text
.ft CW
.nf
.ne \\$1
..
.de Ve \" End verbatim text
.ft R
.fi
..
.\" Set up some character translations and predefined strings.  \*(-- will
.\" give an unbreakable dash, \*(PI will give pi, \*(L" will give a left
.\" double quote, and \*(R" will give a right double quote.  \*(C+ will
.\" give a nicer C++.  Capital omega is used to do unbreakable dashes and
.\" therefore won't be available.  \*(C` and \*(C' expand to `' in nroff,
.\" nothing in troff, for use with C<>.
.tr \(*W-
.ds C+ C\v'-.1v'\h'-1p'\s-2+\h'-1p'+\s0\v'.1v'\h'-1p'
.ie n \{\
.    ds -- \(*W-
.    ds PI pi
.    if (\n(.H=4u)&(1m=24u) .ds -- \(*W\h'-12u'\(*W\h'-12u'-\" diablo 10 pitch
.    if (\n(.H=4u)&(1m=20u) .ds -- \(*W\h'-12u'\(*W\h'-8u'-\"  diablo 12 pitch
.    ds L" ""
.    ds R" ""
.    ds C` ""
.    ds C' ""
'br\}
.el\{\
.    ds -- \|\(em\|
.    ds PI \(*p
.    ds L" ``
.    ds R" ''
.    ds C`
.    ds C'
'br\}
.\"
.\" Escape single quotes in literal strings from groff's Unicode transform.
.ie \n(.g .ds Aq \(aq
.el       .ds Aq '
.\"
.\" If the F register is >0, we'll generate index entries on stderr for
.\" titles (.TH), headers (.SH), subsections (.SS), items (.Ip), and index
.\" entries marked with X<> in POD.  Of course, you'll have to process the
.\" output yourself in some meaningful fashion.
.\"
.\" Avoid warning from groff about undefined register 'F'.
.de IX
..
.nr rF 0
.if \n(.g .if rF .nr rF 1
.if (\n(rF:(\n(.g==0)) \{\
.    if \nF \{\
.        de IX
.        tm Index:\\$1\t\\n%\t"\\$2"
..
.        if !\nF==2 \{\
.            nr % 0
.            nr F 2
.        \}
.    \}
.\}
.rr rF
.\"
.\" Accent mark definitions (@(#)ms.acc 1.5 88/02/08 SMI; from UCB 4.2).
.\" Fear.  Run.  Save yourself.  No user-serviceable parts.
.    \" fudge factors for nroff and troff
.if n \{\
.    ds #H 0
.    ds #V .8m
.    ds #F .3m
.    ds #[ \f1
.    ds #] \fP
.\}
.if t \{\
.    ds #H ((1u-(\\\\n(.fu%2u))*.13m)
.    ds #V .6m
.    ds #F 0
.    ds #[ \&
.    ds #] \&
.\}
.    \" simple accents for nroff and troff
.if n \{\
.    ds ' \&
.    ds ` \&
.    ds ^ \&
.    ds , \&
.    ds ~ ~
.    ds /
.\}
.if t \{\
.    ds ' \\k:\h'-(\\n(.wu*8/10-\*(#H)'\'\h"|\\n:u"
.    ds ` \\k:\h'-(\\n(.wu*8/10-\*(#H)'\`\h'|\\n:u'
.    ds ^ \\k:\h'-(\\n(.wu*10/11-\*(#H)'^\h'|\\n:u'
.    ds , \\k:\h'-(\\n(.wu*8/10)',\h'|\\n:u'
.    ds ~ \\k:\h'-(\\n(.wu-\*(#H-.1m)'~\h'|\\n:u'
.    ds / \\k:\h'-(\\n(.wu*8/10-\*(#H)'\z\(sl\h'|\\n:u'
.\}
.    \" troff and (daisy-wheel) nroff accents
.ds : \\k:\h'-(\\n(.wu*8/10-\*(#H+.1m+\*(#F)'\v'-\*(#V'\z.\h'.2m+\*(#F'.\h'|\\n:u'\v'\*(#V'
.ds 8 \h'\*(#H'\(*b\h'-\*(#H'
.ds o \\k:\h'-(\\n(.wu+\w'\(de'u-\*(#H)/2u'\v'-.3n'\*(#[\z\(de\v'.3n'\h'|\\n:u'\*(#]
.ds d- \h'\*(#H'\(pd\h'-\w'~'u'\v'-.25m'\f2\(hy\fP\v'.25m'\h'-\*(#H'
.ds D- D\\k:\h'-\w'D'u'\v'-.11m'\z\(hy\v'.11m'\h'|\\n:u'
.ds th \*(#[\v'.3m'\s+1I\s-1\v'-.3m'\h'-(\w'I'u*2/3)'\s-1o\s+1\*(#]
.ds Th \*(#[\s+2I\s-2\h'-\w'I'u*3/5'\v'-.3m'o\v'.3m'\*(#]
.ds ae a\h'-(\w'a'u*4/10)'e
.ds Ae A\h'-(\w'A'u*4/10)'E
.    \" corrections for vroff
.if v .ds ~ \\k:\h'-(\\n(.wu*9/10-\*(#H)'\s-2\u~\d\s+2\h'|\\n:u'
.if v .ds ^ \\k:\h'-(\\n(.wu*10/11-\*(#H)'\v'-.4m'^\v'.4m'\h'|\\n:u'
.    \" for low resolution devices (crt and lpr)
.if \n(.H>23 .if \n(.V>19 \
\{\
.    ds : e
.    ds 8 ss
.    ds o a
.    ds d- d\h'-1'\(ga
.    ds D- D\h'-1'\(hy
.    ds th \o'bp'
.    ds Th \o'LP'
.    ds ae ae
.    ds Ae AE
.\}
.rm #[ #] #H #V #F C
.\" ========================================================================
.\"
.IX Title "run-saotrace 1"
.TH run-saotrace 1 "2020-10-06" "2.0.5" "User Contributed Perl Documentation"
.\" For nroff, turn off justification.  Always turn off hyphenation; it makes
.\" way too many mistakes in technical documents.
.if n .ad l
.nh
.SH "NAME"
.Vb 1
\& run\-saotrace \- run SAOTrace commands in a container
.Ve
.SH "SYNOPSIS"
.IX Header "SYNOPSIS"
.Vb 1
\& run\-saotrace [options] <command> [options for <command>]
\&
\& # run a ray trace using docker
\& run\-saotrace \-l docker trace\-nest
\&
\& # look at the ray trace parameters
\& run\-saotrace \-l docker plist trace\-nest
\&
\& # jump into the container and run a bash shell
\& run\-saotrace \-l docker bash
.Ve
.SH "DESCRIPTION"
.IX Header "DESCRIPTION"
\&\fBsaotrace\fR runs the specified command in a container, passing it the
optional arguments.  It can use either \f(CW\*(C`docker\*(C'\fR or \f(CW\*(C`podman\*(C'\fR to
launch the container. If \f(CW\*(C`docker\*(C'\fR is used and the user is not in the
\&\f(CW\*(C`docker\*(C'\fR group, \f(CW\*(C`docker\*(C'\fR will be launched with \f(CW\*(C`sudo\*(C'\fR.
.SH "OPTIONS AND ARGUMENTS"
.IX Header "OPTIONS AND ARGUMENTS"
.IP "\-\-launcher | \-l \fIlauncher\fR" 4
.IX Item "--launcher | -l launcher"
Specify the launcher which will start the container.  It may be one of
\&\f(CW\*(C`docker\*(C'\fR, \f(CW\*(C`podman\*(C'\fR, or \f(CW\*(C`any\*(C'\fR.  If \f(CW\*(C`any\*(C'\fR, \f(CW\*(C`podman\*(C'\fR will be used if
available, otherwise \f(CW\*(C`docker\*(C'\fR.  The \s-1RUN_SAOTRACE_LAUNCHER\s0
environment variable may be used as well.
.Sp
It defaults to \f(CW\*(C`any\*(C'\fR.
.IP "\-\-workdir=\fI\s-1DIR\s0\fR" 4
.IX Item "--workdir=DIR"
cd to \fI\s-1DIR\s0\fR in the container before starting the command. defaults to
the current directory.
.IP "\-\-image=\fI\s-1IMAGE\s0\fR" 4
.IX Item "--image=IMAGE"
Use the specified image to create the container.  It defaults to an appropriate one.
.IP "\-\-help" 4
.IX Item "--help"
output short usage instructions and exit
.IP "\-\-manual" 4
.IX Item "--manual"
output a manual page and exit
.IP "\-\-version" 4
.IX Item "--version"
output version informaiton and exit
.SH "CONTAINER"
.IX Header "CONTAINER"
The official SAOTrace container packages the SAOTrace executables as well as the
official Chandra raytrace configuration databases.  No external resources are required,
except for the \f(CW\*(C`trace\-nest\*(C'\fR parameter file.  This can be created in the work directory
via
.PP
.Vb 1
\& run\-saotrace pset trace\-nest config_dir=\*(Aq${SAOTRACE_DB}/ts_config\*(Aq config_db=\*(Aq${CONFIGDB}\*(Aq
.Ve
.PP
(Mind those single quotes!)
.SS "Container Environment"
.IX Subsection "Container Environment"
The following environment variables are set
.IP "\s-1PFILES\s0" 4
.IX Item "PFILES"
This is set so that parameter files will be written to the working
directory.  It is not possible to cause them to be written elsewhere.
.IP "\s-1SAOTRACE_DB\s0" 4
.IX Item "SAOTRACE_DB"
This is set to point at the root of the Chandra raytrace configuration stored in the container.
.IP "\s-1CONFIGDB\s0" 4
.IX Item "CONFIGDB"
This is set to the official Chandra raytrace configuration, which is currently \f(CW\*(C`orbit\-200809\-01f\-a\*(C'\fR.
.SH "ENVIRONMENT"
.IX Header "ENVIRONMENT"
The following environment variables are available:
.PP
\fI\s-1RUN_SAOTRACE_LAUNCHER\s0\fR
.IX Subsection "RUN_SAOTRACE_LAUNCHER"
.PP
The may be set to \f(CW\*(C`any\*(C'\fR, \f(CW\*(C`docker\*(C'\fR, \f(CW\*(C`podman\*(C'\fR.  This is equivalent to
setting the command line option \f(CW\*(C`\-\-launcher\*(C'\fR.
.SH "COPYRIGHT AND LICENSE"
.IX Header "COPYRIGHT AND LICENSE"
This software is copyright the Smithsonian Astrophysical Observatory
and is released under the \s-1GNU\s0 General Public License.  You may find a
copy at <http://www.fsf.org/copyleft/gpl.html>.
.SH "AUTHOR"
.IX Header "AUTHOR"
Diab Jerius <djerius@cfa.harvard.edu>
# >>> roff

=pod

=head1 NAME

 run-saotrace - run SAOTrace commands in a container

=head1 SYNOPSIS

 run-saotrace [options] <command> [options for <command>]

 # run a ray trace using docker
 run-saotrace -l docker trace-nest

 # look at the ray trace parameters
 run-saotrace -l docker plist trace-nest

 # jump into the container and run a bash shell
 run-saotrace -l docker bash

=head1 DESCRIPTION

B<saotrace> runs the specified command in a container, passing it the
optional arguments.  It can use either C<docker> or C<podman> to
launch the container. If C<docker> is used and the user is not in the
C<docker> group, C<docker> will be launched with C<sudo>.

=head1 OPTIONS AND ARGUMENTS

=over

=item --launcher | -l I<launcher>

Specify the launcher which will start the container.  It may be one of
C<docker>, C<podman>, or C<any>.  If C<any>, C<podman> will be used if
available, otherwise C<docker>.  The L<RUN_SAOTRACE_LAUNCHER>
environment variable may be used as well.

It defaults to C<any>.

=item --workdir=I<DIR>

cd to I<DIR> in the container before starting the command. defaults to
the current directory.

=item --image=I<IMAGE>

Use the specified image to create the container.  It defaults to an appropriate one.

=item --help

output short usage instructions and exit

=item --manual

output a manual page and exit

=item --version

output version informaiton and exit

=back

=head1 CONTAINER

The official SAOTrace container packages the SAOTrace executables as well as the
official Chandra raytrace configuration databases.  No external resources are required,
except for the C<trace-nest> parameter file.  This can be created in the work directory
via

 run-saotrace pset trace-nest config_dir='${SAOTRACE_DB}/ts_config' config_db='${CONFIGDB}'

(Mind those single quotes!)

=head2 Container Environment

The following environment variables are set

=over

=item PFILES

This is set so that parameter files will be written to the working
directory.  It is not possible to cause them to be written elsewhere.

=item SAOTRACE_DB

This is set to point at the root of the Chandra raytrace configuration stored in the container.

=item CONFIGDB

This is set to the official Chandra raytrace configuration, which is currently C<orbit-200809-01f-a>.

=back


=head1 ENVIRONMENT

The following environment variables are available:


=head3 RUN_SAOTRACE_LAUNCHER

The may be set to C<any>, C<docker>, C<podman>.  This is equivalent to
setting the command line option C<--launcher>.

=head1 COPYRIGHT AND LICENSE

This software is copyright the Smithsonian Astrophysical Observatory
and is released under the GNU General Public License.  You may find a
copy at L<http://www.fsf.org/copyleft/gpl.html>.

=head1 AUTHOR

Diab Jerius E<lt>djerius@cfa.harvard.eduE<gt>
